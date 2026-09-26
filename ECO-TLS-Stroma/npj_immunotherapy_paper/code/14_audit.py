"""Style audit: burstiness, AI-tic phrases, transition density, sentence uniformity."""
import re
import numpy as np

md = open("/home/user/npj_immunotherapy_paper/manuscript/manuscript.md").read()
# prose only: drop refs, tables, legends, headers
prose = md.split("## References")[0]
prose = re.sub(r"\|.*", "", prose)
prose = re.sub(r"^#+.*", "", prose, flags=re.M)
prose = re.sub(r"\*\*", "", prose)
prose = re.sub(r"\[[0-9,\-–\s]+\]", "", prose)
prose = re.sub(r">.*", "", prose)
sents = [s.strip() for s in re.split(r"(?<=[.!?])\s+", prose) if len(s.split()) > 3]
lens = np.array([len(s.split()) for s in sents])
print(f"sentences: {len(sents)} | mean len: {lens.mean():.1f} | std: {lens.std():.1f} | CV (burstiness): {lens.std()/lens.mean():.2f}  [human target >0.55]")
print(f"short (<=8w): {(lens<=8).sum()} | long (>=35w): {(lens>=35).sum()}")
paras = [p for p in prose.split("\n\n") if len(p.split()) > 10]
pl = np.array([len(p.split()) for p in paras])
print(f"paragraphs: {len(paras)} | para words mean: {pl.mean():.0f} std: {pl.std():.0f}")

AI_TICS = ["delve", "tapestry", "landscape", "crucial", "vibrant", "notably", "importantly",
           "moreover", "furthermore", "in conclusion", "it is important to note", "it is worth noting",
           "plays a key role", "play a key role", "sheds light", "paves the way", "bustling",
           "in today's", "ever-evolving", "additionally", "consequently", "nevertheless",
           "a wide range of", "a variety of", "boasts", "showcase", "underscore", "underscores",
           "multifaceted", "intricate", "meticulous", "robust", "leverage", "harness",
           "testament to", "rich tapestry", "deep dive", "overall,", "in summary", "firstly",
           "secondly", "thirdly", "on the other hand", "in contrast,", "paradigm"]
low = prose.lower()
hits = {t: len(re.findall(r"\b" + re.escape(t) + r"\b", low)) for t in AI_TICS}
hits = {k: v for k, v in hits.items() if v}
print(f"\nAI-tic hits ({sum(hits.values())} total): {hits}")
emd = prose.count("—") + prose.count(" – ")
print(f"em-dashes: {emd}  [human target: few]")
starts = {}
for s in sents:
    w = s.split()[0].strip("(*\"").lower()
    starts[w] = starts.get(w, 0) + 1
rep = {k: v for k, v in sorted(starts.items(), key=lambda x: -x[1])[:8]}
print(f"top sentence openers: {rep}")
print(f"'we ' count: {len(re.findall(r'\\bwe\\b', low))} | 'our ' count: {len(re.findall(r'\\bour\\b', low))}")
# triplet pattern: X, Y and Z / X, Y, and Z frequency
trip = len(re.findall(r"\w+, \w+ (,|and )", prose))
print(f"serial-list patterns: {trip}")
