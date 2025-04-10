#!/usr/bin/env python

with open("./include/acronyms.source", "r") as f:
    lines = f.readlines()

split_lines = []
for line in lines:
    acr, *long_form = line.split()
    split_lines.append([acr, " ".join(long_form)])

split_lines.sort(key=lambda x: x[0])
with open("./include/acronyms.source", "w") as f:
    for acr, long_form in split_lines:
        f.write(f"{acr:<10} {long_form}\n")