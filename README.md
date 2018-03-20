Valerio Ippolito - INFN Sezione di Roma

To format code based on the `.clang-format` file in this folder:

> make sure clang-format is installed; on MacOSX, do:
> ```/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"```
> ```brew install clang-format```

```clang-format -i *.h *.C```

Please make sure you restore the absence of spaces in the `R__LOAD_LIBRARY` macro calls used by runAsymptoticsCLs.C.

Debug level is defined as
  * 0 = verbose
  * 1 = debug
  * 2 = warning
  * 3 = error
  * 4 = fatal
  * 5 = silent
