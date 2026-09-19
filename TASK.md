# EPISEARCH → modern C++

A two-week, ~15–20 hour project on your own 2018 masters code. The goal is not
to learn C++ from scratch. It is to be able to say "yes" to *do you know C++*
and then have somewhere to point.

**Read this first:** the code in `srcf/` is not modern C++ and never was. It is
C-with-`std::vector`: raw pointers as return values, function-local `static`
buffers, `using namespace std;` at namespace scope, 16-parameter functions,
no tests, no build system worth the name. That is the raw material. Every
stage below turns one part of it into something you would put in a pull
request at a trading firm, and — more to the point — leaves you with a
sentence you can say out loud about *why*.

## Ground rules

1. Branch: `git checkout -b modernise`. One commit per stage, minimum.
2. Keep a `NOTES.md`. Each stage ends with a paragraph: what you changed and
   **why**. The notes are the actual interview prep. Nobody will ask whether it
   compiles; they will ask why you chose that over the alternative.
3. When stuck, give it 45 minutes before you come back to me. Then ask — I'll
   give you the next question, not the answer.
4. Don't skip Stage 6. The write-up is what converts the work into interview
   material.

## Setup (macOS)

**Run every command in this document from the repository root** — the folder
containing `CMakeLists.txt` and `srcf/`. Nothing is ever run from inside
`build/`. If you remember `mkdir build && cd build && cmake ..` from 2018,
presets replace that: `cmake --preset dev` creates `build/dev/` for you and
`cmake --build --preset dev` finds it again by name.

First, extract the scaffold into the repo root if you haven't already, then:

```bash
brew install cmake libomp    # both bottled - seconds, not minutes
cmake --preset dev
cmake --build --preset dev -j
ctest --preset dev
```

**Do not `brew install llvm`.** Apple clang handles everything through Stage 3.
Homebrew's llvm is a ~3 GB formula and if there's no bottle for your macOS it
builds from source, which takes hours. The only thing it buys you is
`<stop_token>` in Stage 4, and there's a workaround for that when you get
there.

Without `libomp` the `#pragma omp` lines are ignored and everything runs
serial. That's fine for Stages 0–3; you'll need threading properly by Stage 4,
and by then you'll have deleted the pragmas anyway.

Presets: `dev` (debug), `asan` (address+UB sanitizer), `tsan` (thread
sanitizer), `strict` (warnings are errors), `release` (optimised, for
benchmarking).

---

## Stage 0 — Baseline · 30 min

```bash
./tools/make_dataset.py --nsamples 300 --nsnps 40 --out data.txt --seed 7
# note the planted triple it prints
./build/dev/episearch -inputfile data.txt -outputfile out.txt \
    -nsamples 300 -nvars 41 -printall n -assoclevel 0.05 -numThread 1
```

The generator plants one genuine three-way interaction: three SNP columns that
are individually uninformative about the phenotype but jointly determine it.
That triple is your ground truth for the rest of the project — a correct search
ranks it far above everything else.

**Deliverable:** you will observe something. Write down exactly what, with the
command line that produced it.

---

## Stage 1 — The crash · 2–3 h

The binary segfaults. Not on some exotic edge case — on the example data that
ships in `Examples/fullSearch/`, and on anything the generator produces, on both
the full-search and the indexed code path. It presumably worked on your 2018
toolchain.

- Reproduce it reliably.
- Find it. **The `asan` preset will not catch this one.** Working out *why a
  sanitizer is blind to it* is half the lesson — that fact narrows the space of
  possible causes a lot.
- Fix it minimally. The fix is small. If you find yourself rewriting the
  function, you haven't found it yet. Do not refactor anything at this stage.

**Acceptance:** `episearch` runs to completion on the 40-SNP dataset, and the
planted triple has by far the highest `IGStrict` of any triple in the output.

**In `NOTES.md`:** what the defect was; why it produces a crash rather than a
quietly wrong number; why it might have appeared to work eight years ago; why
the sanitizer can't see it.

> This stage is the exercise your friend got. Time yourself. Being able to say
> "I found a segfault in my own eight-year-old code and here's the root cause"
> is worth more than any amount of syntax recall.

---

## Stage 2 — Make the compiler your ally · 2–3 h

The build currently emits several dozen warnings under the project's warning
set (see `CMakeLists.txt` — `-Wconversion`, `-Wshadow`, `-Wold-style-cast` and
friends). Get to zero, then switch on the `strict` preset and keep it on for
the rest of the project.

Also in scope:

- `using namespace std;` at namespace scope. Remove every one. Be able to say
  why this is a defect and not a style preference.
- Include-what-you-use. Several files include `<map>` and `<ctime>` for nothing
  and rely on transitive includes for what they do use.
- Header guards → `#pragma once`, or fix the ones that are wrong.
- `const` everywhere it belongs, starting with parameters.
- Signed/unsigned comparisons: **fix the types**, don't paper over them with a
  cast. A cast that silences a warning without changing behaviour is a lie.
- Run `clang-tidy` (config is in `.clang-tidy`, deliberately narrow).

**Acceptance:** `cmake --preset strict && cmake --build --preset strict` is
clean; tests still pass; output on the Stage 1 dataset is byte-identical.

---

## Stage 3 — Ownership and interfaces · 5–6 h · the big one

Rewrite the data and entropy layer into something you'd defend in review.

- `calculateMeasures` returns `double*`. Decide what it *should* return and
  why. While you're there: five unnamed doubles addressed as
  `*(measureArray + 3)` is not a return type, it's a riddle. Name them.
- More than one function hands the caller back storage whose lifetime the
  signature does not describe. Find **all** of them. Replace each with
  something where the ownership is obvious from the declaration alone.
- The dataset is a flat `vector<int>` with column-major indexing spelled out by
  hand at every use site (`d[s + i * nsamples]`). Wrap it. The wrapper should
  make an out-of-bounds access hard to write by accident and should cost
  nothing at runtime in a release build. Be ready to justify that second claim.
- Stringly- and int-typed parameters: `int c` selecting an estimator, `bool
  printall`, `"y"`/`"n"` strings from the command line. Give them types.
- Functions with 12–16 positional parameters, several of them adjacent
  same-type `int`s that are trivially swappable at a call site. Fix that.
- Add tests as you go — the list at the bottom of `tests/test_entropy.cpp` is
  your starting point.

**Constraints for this stage:**

- no `new`, no `delete`, no owning raw pointer anywhere
- no function-local `static` holding mutable state
- anything not mutated is taken by `const&`, or by value where that's cheaper
- pass views, not `(pointer, length)` pairs

**Acceptance:** identical numeric output to Stage 2 on the same dataset;
`strict` still clean; test count meaningfully higher.

**Be ready to talk about:** rule of zero; value vs reference semantics; when a
view is the *wrong* choice (lifetime); RVO and why returning by value isn't
slow; `std::optional` / `std::expected` for the file readers, which currently
fail silently on a missing file.

---

## Stage 4 — Concurrency without OpenMP · 3–4 h

Run the same input with `-numThread 1` and `-numThread 4` and compare
`means.txt`. The numbers differ. **Work out why before you change anything** —
there is more than one reason, and they are not the same kind of problem as
each other.

Then delete the OpenMP pragmas and do it in standard C++: `std::jthread`,
per-thread accumulation, an explicit merge step. (If Apple clang's standard
library doesn't give you `<stop_token>` / `std::jthread`, either build with the
Homebrew `llvm` toolchain or write the RAII join-on-destruction wrapper around
`std::thread` yourself — that's a worthwhile ten minutes either way, and it's a
question they ask.)

**Requirements:**

- identical output between 1 and N threads — you define the ordering that makes
  that true, and say so in the notes
- clean under the `tsan` preset
- no mutex in the inner loop

**Be ready to talk about:** data race vs race condition; false sharing and how
your accumulator layout avoids it; why floating-point summation order changes
the result and what you did about it; why you didn't reach for
`std::atomic<double>`; where `std::execution::par` would and wouldn't fit.

---

## Stage 5 — Make it fast · 2–3 h

Benchmark first:

```bash
cmake --preset release && cmake --build --preset release -j
./tools/bench.sh build/release/episearch data.txt 300 41
```

Record the numbers. A speedup you can't quantify is a speedup you can't defend.

The inner loop estimates entropies over a handful of columns, and it does that
tens of millions of times. Profile it — Instruments on macOS, `perf` on Linux.
One data-structure decision dominates the runtime. Find it by measurement, not
by guessing, then change it and measure again.

**Target:** 10× single-threaded on the same dataset, tests still green, output
unchanged.

**Be ready to talk about:** your before/after numbers; what the profile
actually showed; allocation in hot loops; cache behaviour; why the complexity
class was never the problem here.

---

## Stage 6 — The write-up · 1 h

Rewrite `README.md` for the current state of the code. Open a PR from
`modernise` against your own `main` and write the description as if a reviewer
at Citadel Securities were reading it: what changed, why, what you deliberately
did not do.

Then make sure you can answer, out loud, without notes:

- rule of zero / three / five — and which one this codebase now follows
- move semantics: what makes a move cheap, when the compiler elides it entirely
- the cost of virtual dispatch, and why there isn't any in your version
- function-local `static`s: initialisation, thread safety, and why they're a
  trap in a parallel program
- `unique_ptr` vs `shared_ptr` vs neither — for this codebase the answer is
  *neither*, so be able to say why that's a strength
- what undefined behaviour actually licenses the optimiser to do

---

## And the answer to the question itself

Last year: "no." This year, when they ask:

> "I wrote and shipped a C++ research tool for my masters — an information-
> theoretic search over genetic data, parallelised, it ran ~50× faster than the
> R equivalent. I came back to it recently and modernised it to C++20: found and
> fixed a latent segfault that a newer toolchain exposed, replaced the OpenMP
> layer with `std::jthread`, took the raw-pointer and static-buffer interfaces
> out, put it under CMake with sanitizer builds and a test suite, and got a 10×
> on the hot path. I'd call myself competent in C++ rather than expert — I'm
> strongest in Python and Java day to day — but I can read it, debug it, and
> write correct modern C++."

That's a true sentence you will be able to say in two weeks. It is a much
better answer than "no", and — unlike "yes" on its own — every clause in it is
something you can be drilled on and survive.
