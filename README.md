# ELOSA: Exhaustive Lesion Overlap Space Analysis 



ELOSA finds every **group**. A group is a set of patients whose lesions all overlap somewhere, such that no other patient's lesion covers the whole of that overlap region. Each group corresponds to exactly one distinct region that the dataset can tell apart. So the number of groups is the ceiling on what the lesion data can express.

This problem is the same as **closed-itemset mining** (also called formal concept analysis):

* **Transactions:** voxels with identical patterns of lesioned patients, merged into weighted "atoms".
* **Items:** patients.
* **Groups:** the closed itemsets.

The fastest known methods for this problem are the LCM family (Uno et al., 2004). `src/elosa.cpp` is a multithreaded C++17 implementation of that approach. MATLAB functions around it prepare the data, run the program, read the results, and carry out the downstream critical-lesion-site search. The binary lesion data that drive the programme cannot be shared, but if you run it with your own data,  it should be represented as a patients x voxels 2D matrix. Optionally, you can include a vector of impairment of labels, which allows the programme to calculate how predictive is damage in each group, for that impairment (using PPV).

>
> Please run `runtests('ElosaTest')` first.

## Layout

|Path|Contents|
|-|-|
|`src/elosa.cpp`|The enumerator: one file with no dependencies.|
|`matlab/`|Runs the program from MATLAB, reads its output, and holds the analysis functions (see below).|
|`benchmark/`|`elosa\_benchmark.m` times the new program against the original `FindRecursive` on your own data. It includes a copy of the original `ExhaustiveRegionFinder` with only its class name fixed.|
|`tests/`|`ElosaTest.m` (MATLAB) and the Python tests and benchmarks used during development.|
|`archive/`|Original files that weren't ported (see `archive/README.md`).|
|`CMakeLists.txt`|Optional CMake build.|

## Building

The program must be compiled once on each machine. The build uses `-march=native`, so the binary is tuned to the CPU it's built on.

* **From MATLAB with MinGW:** install the *MATLAB Support for MinGW-w64 C/C++ Compiler* add-on, then run:

```matlab
  addpath matlab
  elosa\_compile          % builds bin\\elosa.exe
  ```

  `elosa\_compile` tries the fastest settings first, then safer ones:

  * a generic CPU instead of `-march=native`;
  * Windows' own threads, for MinGW builds without `std::thread` (the "win32" thread model);
  * a non-static build.

  Everything the compiler prints goes to `bin\\compile\_log.txt`. If all the attempts fail, the error message shows that output.

* **Visual Studio:** open a *Developer Command Prompt* in this folder and run:

```
  cl /O2 /std:c++17 /EHsc /Fe:bin\\elosa.exe src\\elosa.cpp
  ```

* **Linux or macOS:**

```
  g++ -O3 -march=native -std=c++17 -pthread -o bin/elosa src/elosa.cpp
  ```

  Or use `cmake -S . -B build \&\& cmake --build build`.

## Quick start

```matlab
addpath matlab
S = elosa\_select\_patients(D, {'L1==English', 'RightHanded\_Only', 'FirstScanOnly'});
Impaired = logical(S.Labels(:, 3));

R = elosa\_run(S.Lesions, 'Use', Impaired);                 % count only: fast, nothing stored
R.count, R.by\_size

R = elosa\_run(S.Lesions, 'Use', Impaired, 'Covering', true);   % every group + who covers it
Solution = elosa\_greedy\_cls(R.Covering, Impaired, 'MinPPV', 0.9, 'RegionVoxels', R.RegionVoxels);
Regions  = elosa\_regions(R.Groups(Solution.Group, :), S.Lesions);   % voxel masks
```

`matlab/elosa\_example.m` walks through a complete analysis, including a train/test split and a growth curve. `D` can be one of three things:

* a `CLS` object;
* a `PLORAS\_parallel` object;
* the struct from the Confidence project's `load\_ploras` (it uses the same field names).

|Function|Replaces|Purpose|
|-|-|-|
|`elosa\_run`|`ExhaustiveRegionFinder`, the `ELOSA\_v\*` drivers, `RunRecursiveParallel\*`, `GrowGroup\*`, `Mex\*Find`|Counts or lists every group. It can also estimate the count.|
|`elosa\_select\_patients`|`GetGoodData`|Applies the inclusion criteria.|
|`elosa\_regions`|`prod(Lesions(group,:))` loops|Region masks for groups.|
|`elosa\_overlaps`|`OverlapLesions`|Fractional overlaps. Only needed for thresholds below 1.|
|`elosa\_greedy\_cls`|`GetCLS`, `FindGreedyCLS`, `GetUnexplained`|Greedy critical-lesion-site search, with optional train/test evaluation.|
|`elosa\_growth\_curve`|`MakeEOLSAGrowthCurve`|Group counts for random subsets of increasing size.|
|`elosa\_downsample`|`ReduceResolution`, `ChangeLesionResolutios`|Coarser voxels in pure MATLAB, with no SPM needed.|
|`elosa\_write\_input`, `elosa\_read\_groups`, `elosa\_compile`, `elosa\_exe`, `elosa\_system`, `elosa\_mingw\_bin`, `elosa\_progress\_display`, `elosa\_stop`|–|Plumbing.|

`Use` marks the patients who may form groups, usually the impaired ones. The other patients don't define groups, but the program counts them whenever their lesion covers a region. So each group comes with:

* `NumOther`: how many non-grouped patients cover the region, which gives the PPV at threshold 1;
* `Covering`: which patients those are (optional).

This replaces the patients × groups `Overlaps` matrix, which was the other big memory and time cost.

## Progress while it runs

`elosa\_run` shows a status line in the Command Window for any run longer than a second. It updates in place every half second:

```
ELOSA: 1,234,567 groups so far, 2m 10s elapsed, roughly 34% done, 3m 10s to 6m 20s left
```

* **Exact:** the number of groups found so far and the time elapsed.
* **Rough:** the percentage and the time left. The search tree's branches differ hugely in size, so the program can only judge what's left from the branches it has finished. On test data its guess of the total time was 0.9–2 times the real time, so the time left is shown as a range, from half its guess to its full guess.
* **Other phases:** reading the input and finding the voxel patterns are shown too, and `'Estimate'` runs show how many of their random walks are done.

You can change the display:

* `'Progress', 'bar'` shows a progress bar with a Cancel button.
* `'Progress', 'none'` turns the display off. `elosa\_growth\_curve` uses this by default and prints one line per run instead.

**Stopping a run:**

* **Ctrl+C in MATLAB**, or Cancel on the progress bar, stops the program within about half a second. Earlier versions left it running in the background.
* **If MATLAB is closed or crashes**, the program notices within a second and stops by itself (`--parent-pid`).
* **`elosa\_stop`** ends any copies still running from before this change. So does `taskkill /F /IM elosa.exe` in a Command Prompt.
* **An out-of-date build** is refused: `elosa\_run` stops if `bin\\elosa.exe` is older than `src\\elosa.cpp` and asks you to run `elosa\_compile`. When run from the command line, `elosa` prints the same information every 30 seconds, and `--progress FILE` rewrites FILE every second.

### The space grows exponentially

At 150–200 patients with realistic overlap, it can reach billions of groups, and no implementation will enumerate that quickly. Here are the tools for that situation:

* **Count mode** (the default) is the quickest option, and is all you need for the ceiling.
* **`'Estimate', N`** estimates the count from N random walks (Knuth's estimator), usually in seconds. It is unbiased in principle and was accurate on small spaces (e.g. 3,899 ± 106 for a true 3,875). On large spaces it underestimates in practice, because the tree is very uneven: it gave 1.6–2.1 million for a true 3.3 million. **Treat the estimate as an order of magnitude.**
* **`elosa\_growth\_curve`** gives exact counts on subsets of increasing size, so you can extrapolate.
* **`'MinVoxels'`** counts only regions of a minimum size, which prunes the search. The old code computed `MinSize = floor(800/prod(VoxSize))` but never used it.
* **`'MaxPatients'`** caps the group size, which also prunes the search.
* **`elosa\_downsample`** gives coarser voxels, so there are fewer distinct patterns.

## Tests

* **`tests/test\_vs\_brute.py`** runs thousands of random datasets. Each dataset checks every group, region size, covering count and covering list against a brute-force implementation of the definition. The runs cover 1–4 threads, all three orderings, `--min-voxels` and `--max-patients`.
* **`tests/test\_large.py`** covers 65–200 patients, including multi-word bitsets, non-grouped patients and thread hand-offs.
* **Sanitizers:** both suites also passed with AddressSanitizer/UBSan and ThreadSanitizer builds.
* **`tests/original\_port.py`** is the NumPy port of the original search. The missing and bogus groups above come from it.
* **`tests/ElosaTest.m`** is the MATLAB test suite. 

