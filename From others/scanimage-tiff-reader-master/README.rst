.. image:: https://gitlab.com/vidriotech/scanimage-tiff-reader/badges/master/pipeline.svg
   :target: https://gitlab.com/vidriotech/scanimage-tiff-reader/commits/master
   :alt: Pipeline status

.. image:: https://gitlab.com/vidriotech/scanimage-tiff-reader/badges/master/coverage.svg
   :target: https://gitlab.com/vidriotech/scanimage-tiff-reader/commits/master
   :alt: Coverage

About
=====

For more information, see the documentation_.

The ScanImageTiffReader is a command line interface and C library that
reads data from Tiff_ and BigTiff_ files recorded using ScanImage_.  It was
written with performance in mind and provides access  to ScanImage-specific
metadata. This library should actually work with most tiff files, but as of
now we don't support compressed or tiled data. 

Cross-platform language bindings are available for Python_ (2.7.14 and 3.6), Matlab_
(2016b+), or Julia_ (1.0).  Those are available from their respective package
managers.  See the links for more information.

Both ScanImage_ and this reader are products of `Vidrio Technologies`_.  If you
have questions or need support feel free to email_.

.. _documentation: https://vidriotech.gitlab.io/scanimage-tiff-reader/
.. _Tiff: https://en.wikipedia.org/wiki/Tagged_Image_File_Format
.. _BigTiff: http://bigtiff.org/
.. _ScanImage: http://scanimage.org
.. _scanimage.org: http://scanimage.org
.. _Matlab: https://vidriotech.gitlab.io/scanimagetiffreader-matlab
.. _Julia: https://vidriotech.gitlab.io/scanimagetiffreader-julia
.. _Python: https://vidriotech.gitlab.io/scanimagetiffreader-python
.. _`Vidrio Technologies`: http://vidriotechnologies.com/
.. _email: support@vidriotech.com

Examples
========

Shell
`````

.. code-block:: bash

    ScanImageTiffReader image data my.tif > my.raw

C/C++
`````

.. code-block:: c

    #include <tiff.reader.api.h>
    const size_t bytes_per_pixel[]={1,2,4,8,1,2,4,8,4,8};
    ScanImageTiffReader reader=ScanImageTiffReader_Open("my.tif");
    struct nd shape=ScanImageTiffReader_GetShape(reader);
    size_t nbytes=bytes_per_pixel[shape->type]*shape.strides[shape.ndim];
    void *data=malloc(nbytes);
    ScanImageTiffReader_GetData(reader,data,nbytes);
    ScanImageTiffReader_Close(reader);

Building
========

The source code is tracked using git_, which also provides a capability for
managing links to other modules.  For building the documentation and running
tests it may be necessary to download these modules.  To do so, open a
console, navigate to the root source directory, and run::

    git submodule update --init

CMake_ is used to configure the build.  Starting in the source directory,
you use it something like this::

    mkdir build
    cd build
    cmake ..

At that point CMake will generate a Makefile or a Visual Studio project or
something else.  It supports the ability to generate project files for many
kinds of build systems.  CMake also has a GUI_ to help configure the
project.

Further, we use several features of cmake to automate testing (CTest_) and
packaging (CPack_).

Testing
```````

To run the tests after building, from the build directory, run::

    ctest

The tests use data in the ``data`` directory of the repository.  If the data
directory isn't there, it may be necessary to run ::

    git submodule update --init

Installation and Packaging
``````````````````````````

CMake_ will generate and install target (e.g. ``make install``).  To control
where files are installed, set the CMAKE_INSTALL_PREFIX_ during configuration.

CPack_ can be used to bundle the installed files into an installer or a
zip-file.  Several options are available and they can be selected by setting
the appropriate 'CPACK_*' variables during configuration.

.. _git: https://git-scm.com/
.. _CMake: https://cmake.org/
.. _CPack: https://cmake.org/cmake/help/v3.5/manual/cpack.1.html
.. _CTest: https://cmake.org/cmake/help/v3.5/manual/ctest.1.html
.. _CMAKE_INSTALL_PREFIX: https://cmake.org/cmake/help/v3.0/variable/CMAKE_INSTALL_PREFIX.html
.. _GUI: https://cmake.org/runningcmake/


Documentation
`````````````

Building the documentation requires Python, sphinx_, and doxygen_.  To
set these up:

1. Install python.  Make sure it's on your PATH.

   On windows, I recommend using Anaconda_ installs python together with a
   number of useful packages.

2. Install doxygen.
3. Install sphinx., from the command line, run::

    pip install sphinx

or::

    pip install -r doc/requirements.txt

.. _sphinx: http://www.sphinx-doc.org/en/stable/
.. _doxygen: http://www.stack.nl/~dimitri/doxygen/
.. _Anaconda: https://www.continuum.io/downloads

Finally, to build the documentation, navigate to the source directory
using the command line, and then run::

    git submodule update --init
    cd doc
    doxygen -u
    make html

If everything goes smoothly, ``index.html`` (the langing page) will be in
``doc/build/html/``.  On windows, the ``make`` command is actually
``make.bat``.
