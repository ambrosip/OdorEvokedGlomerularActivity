Download
~~~~~~~~

Pre-built binaries are available for 64-bit Windows, OS X, or Linux.

=========== ========== ==========
   Target    Version    Download
----------- ---------- ----------
Windows x64    1.4.1     `Download <win-recent_>`_ 
OS X/Linux     1.4.1     We want to hear from you at support@mbfbioscience.com!
=========== ========== ========== 

.. _here: http://scanimage.vidriotechnologies.com/display/SIH/Tools
.. _win-recent: https://gitlab.com/api/v4/projects/922490/jobs/artifacts/master/raw/ScanImageTiffReader-1.4.1-win64.zip?job=build_windows
.. _win-130: https://gitlab.com/vidriotech/scanimage-tiff-reader/-/jobs/artifacts/1.3/download?job=build_windows
.. _osx-130: https://gitlab.com/vidriotech/scanimage-tiff-reader/-/jobs/artifacts/1.3/download?job=build_osx
.. _nix-130: https://gitlab.com/vidriotech/scanimage-tiff-reader/-/jobs/artifacts/1.3/download?job=build_nix

About
~~~~~

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
~~~~~~~~

Each example shows how to extract raw image data from the file.

:doc:`command-line-tool`
------------------------

.. code-block:: bash

    ScanImageTiffReader image data my.tif > my.raw

:doc:`c`
--------

.. code-block:: c

    #include <tiff.reader.api.h>
    const size_t bytes_per_pixel[]={1,2,4,8,1,2,4,8,4,8};
    ScanImageTiffReader reader=ScanImageTiffReader_Open("my.tif");
    struct nd shape=ScanImageTiffReader_GetShape(reader);
    size_t nbytes=bytes_per_pixel[shape->type]*shape.strides[shape.ndim];
    void *data=malloc(nbytes);
    ScanImageTiffReader_GetData(reader,data,nbytes);
    ScanImageTiffReader_Close(reader);
