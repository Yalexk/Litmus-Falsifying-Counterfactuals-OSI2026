"""LLM stage: counterfactual sets -> structured, executable hypotheses.

Uses ANTHROPIC_API_KEY if set, otherwise the `claude` CLI (Claude Code login,
no key needed). The model is told to emit JSON
only. Everything it emits is treated as untrusted input and run through
rules.validate() before it touches the data.
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile

import anthropic

from .rules import GROUP_SMARTS, Rule, validate

MODEL = os.environ.get("CF_MODEL", "claude-sonnet-4-6")

# Bump PROMPT_VERSION with every SYSTEM change; it is recorded in each run's
# manifest. The prompt may describe HOW rules are scored, never WHICH chemistry
# to propose -- whoever edits it has seen test results.
PROMPT_VERSION = "v2"

SYSTEM = """You are generating falsifiable structure-property hypotheses from \
counterfactual explanations of a molecular property model.

You will see a base molecule, its predicted value, and nearby molecules whose \
predicted value differs. Propose general rules that would explain the pattern.

How your rules will be scored, on measured data you will not see:
1. Raw test: do molecules matching the condition differ from the rest in the \
claimed direction, by more than experimental error (~0.6 log units)?
2. Matched-lipophilicity test: the same comparison made only between \
molecules of similar MolLogP. Lipophilicity alone explains most of this \
property, so a rule that is only a proxy for lipophilicity (or for things \
that track it closely, such as overall size or overall polarity) passes the \
raw test and then shows no effect at matched MolLogP. Such rules count as \
restating something already known.
Aim for rules that would still hold between two molecules of the same \
MolLogP. A rule whose condition is nearly equivalent to a MolLogP threshold \
cannot be tested this way and earns nothing.
A condition that matches very few molecules cannot be tested at all; aim for \
conditions that a reasonable fraction of small organic molecules satisfy.

Output JSON only. No prose, no markdown fences. Schema:
{"hypotheses": [
  {"text": "<one sentence claim a chemist could argue with>",
   "direction": "higher" | "lower",
   "condition": {"kind": "descriptor", "name": "<descriptor>", "op": ">"|"<", "value": <number>}
             OR {"kind": "smarts", "pattern": "<valid SMARTS>"}}
]}

Rules:
- "direction" means: molecules MATCHING the condition have higher/lower \
property values than molecules that do not.
- Prefer descriptor conditions. Use SMARTS only for a specific functional group.
- Propose general rules, not restatements of one molecule pair.
- Do not hedge. A vague hypothesis cannot be falsified and is worthless here.
"""

DESCRIPTORS_DOC = """Available descriptors: MolWt, MolLogP, TPSA, NumHDonors, \
NumHAcceptors, NumRotatableBonds, NumAromaticRings, FractionCSP3, RingCount.
Known-good SMARTS for common groups: """ + json.dumps(GROUP_SMARTS)


def _client():
    return anthropic.Anthropic()


def _ask_cli(prompt: str) -> str:
    """No API key: use the Claude Code CLI. It runs from an empty temp dir
    with no tools at all, so the model sees only this prompt -- it cannot
    read the test split, this repo, or its CLAUDE.md, or browse to the
    public ESOL file. `--tools ""` only removes built-in tools; the strict
    empty MCP config is what removes the user's MCP servers (checked via the
    stream-json init event: tools == [])."""
    with tempfile.TemporaryDirectory() as tmp:
        r = subprocess.run(
            ["claude", "-p", "--tools", "", "--strict-mcp-config",
             "--mcp-config", '{"mcpServers":{}}', "--system-prompt", SYSTEM,
             "--model", MODEL, "--no-session-persistence", "--output-format", "text"],
            input=prompt, capture_output=True, text=True, cwd=tmp, timeout=600,
        )
    if r.returncode != 0:
        raise RuntimeError(f"claude CLI failed: {r.stderr.strip()[:500]}")
    return r.stdout


def _ask(messages, max_tokens=2000) -> dict:
    if os.environ.get("ANTHROPIC_API_KEY"):
        r = _client().messages.create(
            model=MODEL, max_tokens=max_tokens, system=SYSTEM, messages=messages
        )
        txt = "".join(b.text for b in r.content if b.type == "text")
    else:
        txt = _ask_cli(messages[-1]["content"])
    txt = txt.strip().removeprefix("```json").removeprefix("```").removesuffix("```")
    return json.loads(txt)


def propose(cf_records: list[dict], property_name: str, n: int = 8) -> list[Rule]:
    """One call over a batch of counterfactual sets -> candidate rules."""
    prompt = (
        f"Property being modelled: {property_name}\n\n{DESCRIPTORS_DOC}\n\n"
        f"Counterfactual evidence:\n{json.dumps(cf_records, indent=1)[:40000]}\n\n"
        f"Propose at most {n} hypotheses."
    )
    data = _ask([{"role": "user", "content": prompt}])
    rules = []
    for h in data.get("hypotheses", []):
        r = Rule(
            text=h.get("text", ""),
            condition=h.get("condition", {}),
            direction=h.get("direction", ""),
            provenance={"stage": "propose"},
        )
        err = validate(r)
        if err:
            print(f"  rejected malformed rule: {err}")
            continue
        rules.append(r)
    return rules


def revise(rule: Rule, result: dict, property_name: str) -> Rule | None:
    """Feed a failed rule back with its evidence. One retry, not a loop --
    unbounded revision just overfits to the validation split."""
    prompt = (
        f"Property: {property_name}\n"
        f"This hypothesis failed held-out testing:\n{json.dumps(result, indent=1)}\n\n"
        f"{DESCRIPTORS_DOC}\n\n"
        "Propose exactly one revised hypothesis. If the original claim is simply "
        "wrong, say so by proposing the opposite direction rather than narrowing "
        "the condition until it fits."
    )
    data = _ask([{"role": "user", "content": prompt}])
    hs = data.get("hypotheses", [])
    if not hs:
        return None
    h = hs[0]
    r = Rule(h.get("text", ""), h.get("condition", {}), h.get("direction", ""),
             {"stage": "revise", "parent": rule.text})
    return None if validate(r) else r
