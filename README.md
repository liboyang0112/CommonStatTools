# Exotics statistics scripts

Valerio Ippolito - INFN Sezione di Roma

based on the work of many people over a decade!

## Introduction
This is a collection of common tools to perform the statistical analysis of a RooWorkspace.
The goal of this repository is not to provide a common framework, but rather to provide
an unique place to find whatever might be needed to go from a workspace to a paper.

These tools try to be as general as possible, and simple to use - whenever this is compatible with
the way they were originally designed, of course... Doxygen has been added and should be regarded
as work in progress: there is much room for improvement, feel free to help!

For support and announcements, join the [atlas-phys-stat-commonstattools](https://e-groups.cern.ch/e-groups/Egroup.do?egroupName=atlas-phys-stat-commonstattools) CERN e-group.

## Quick start
On lxplus:
```
setupATLAS
lsetup git
lsetup "root 6.14.04-x86_64-slc6-gcc62-opt"
lsetup cmake
git clone --recurse-submodules https://:@gitlab.cern.ch:8443/atlas-phys/exot/CommonStatTools.git
cd CommonStatTools
mkdir build
cd build 
cmake ..
make
```

After the first checkout, please make sure you have the latest version with `git pull`. Currently, only the `master` branch is defined.

A testing file is provided, `test.sh`. You can edit it to point to your favourite RooWorkspace and run it.
If your workspace was not generated via HistFactory, some of the tools won't of course run.

To include this package in your macros, you can also simply do something like:
```
#include "CommonStatTools/Minimization.h"

R__LOAD_LIBRARY(CommonStatTools/build/libExoStats.so)
//R__LOAD_LIBRARY(CommonStatTools/build/libExoStats.dylib) // on macOS
// R__LOAD_LIBRARY(CommonStatTools/Minimization.C+) // uncomment this if you did not run cmake+make for some reason

void test() {

RooWorkspace w;
w.factory("Gaussian::pdf(x[-100,100],mean[-10,10],sigma[0.1,10])");
auto data = w.pdf("pdf")->generate(*w.var("x"),100000);
EXOSTATS::fitPdf(w.pdf("pdf"), data);
}
```

## For users

This repository collects code to:
  * calculate limits using profile likelihood, CLs and asymptotic formulae: `runAsymptoticsCLs.C`
  * calculate limits using any test statistics: `StandardHypoTestInvDemo.C`
  * calculate local p0 using profile likelihood and asymptotic formulae: `runSig.C`
  * calculate local p0 using any test statistics: `StandardHypoTestDemo.C`
  * calculate global p0: `getGlobalP0.C`
  * perform cross-checks: `runFitCrossCheckForLimits.C`
  * plot fit covariance matrices: `getCorrMatrix.C`
  * rank nuisance parameters: `StatisticsTools/submit_pulls.py` (batch) and `StatisticsTools/pulls.exe` (local)
  * HistFactory-only: extract histograms from workspace files: `obtainHistosFromWS.py`
  * HistFactory-only: compare histograms between two different workspace files: `compareHistos.py`
  * HistFactory-only: print pre- and post-fit yield tables: `getHFtables.C`
  * HistFactory-only: evaluate the impact of nuisance parameters over a yield (signal or background) in a channel, possibly after fitting different channels: `getHFtables.C`

In a few cases, classes were also made available in order for developers to be able to fully configure and customize these tools: see the dedicated section for more details.

Nuisance parameter ranking relies on a nice tool from Stefan Gadatsch, which is standalone and is included in this repository as a git [submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules). In order to run it, please refer to its `README.md` file, but it's something like
```
cd StatisticsTools
python submit_pulls.py --help # for interactive running
./pulls.exe -h # for local running
```

(Beware in case you have multiple POIs - most standard code, like `runAsymptoticsCLs.C`, works only with the first one!)

## For developers

Fit helper functions are provided, which should be easy to integrate into other frameworks:
  * `Minimization.C` contains functions to minimize functions smartly, mostly taken from Aaron Armbuster's original version of the `runAsymptoticsCLs.C` macro
  * `AsimovDataMaking.C` contains functions to prepare Asimov datasets (also from Aaron's macro)
  * `AsymptoticsCLsRunner.C` implements a class which wraps Aaron's original macro, so that its features can be more easily used by other frameworks

Documenation is based on Doxygen; [ROOT's conventions](https://root.cern.ch/formatting-comments-doxygen) are followed.

Code contained in files whose name starts with a capital letter are tools conceived to be potentially used by / included in other frameworks; otherwise, they are macros or runners (which often rely on the former tools).
 
Almost all tools include a `debugLevel` flag/option, which is however not always implemented fully yet. The debug level is defined as
  * 0 = verbose
  * 1 = debug
  * 2 = warning
  * 3 = error
  * 4 = fatal
  * 5 = silent

To format code please use `clang-format`, which is based on the `.clang-format` file in this folder:

> make sure clang-format is installed; on MacOSX, do:
> ```/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"; ```brew install clang-format```

```clang-format -i *.h *.C```

When you format the code, please make sure you restore the absence of spaces in the `R__LOAD_LIBRARY` macro calls used by macros or runners (e.g. `runAsymptoticsCLs.C`).
