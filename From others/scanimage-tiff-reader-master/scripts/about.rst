Scripts
=======

I use some build scripts to automate building for continuous integration
and sometimes for convenience.

Usually these are designed to run from the root source directory like this::

    powershell -File scripts/build.ps1

=================== ==============================================================
Script              Description
=================== ==============================================================
build.ps1           Builds on Windows.  Assumes Visual Studio 14.0 is installed.
build_init.ps1      Used by ``continuous.ps1`` to do some setup.
build_docs.ps1      Builds documentation on windows.
continuous.ps1      Watches the ``src-cpp`` folder for changes and rebuilds by
                    calling ``build.ps1``.
continuous_docs.ps1 Watches the ``src-cpp`` folder for changes and rebuilds by 
                    calling ``build_docs.ps1``.
=================== ==============================================================
