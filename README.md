# Sample usage
Generate empty `build/` and `obj/` directories:
```bash
make gen
```
Compile:
```bash
make
```
(Optional) Compile and generate `compile_commands.json` with `bear`:
```bash
bear -- make
```
Run:
```bash
make run
```
Test (basic util and gf2 functions):
```bash
make test
```
