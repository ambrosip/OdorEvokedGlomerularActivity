Command Line Interface
======================

.. _ScanImageTiffReader:

ScanImageTiffReader
~~~~~~~~~~~~~~~~~~~

The command line tool is designed to make inspecting files simple.  It also
allows extraction of the image data into a raw file.

.. table:: Basic Commands

    ============================================  ===================================
    Command                                       Description
    ============================================  ===================================
    :option:`ScanImageTiffReader descriptions`    Print all image descriptions.
              with ``descriptions --frame <id>``  Print one image description.
    :option:`ScanImageTiffReader metadata`        Print the :ref:`ScanImage metadata <format>`.
    :option:`ScanImageTiffReader data`            Print the raw image data as a byte stream.
    :option:`ScanImageTiffReader help`            Print detailed help for a
                                                  sub-command.
    ============================================  ===================================

Usage
-----

From the command line::

    ScanImageTiffReader

will output a basic help message::

    Usage: ScanImageTiffReader <command> [<args>]

    These are the available commands:
       descriptions     Extract the contents of a the image description tag(s).
                        If you do not specify which frame, then all the image
                        descriptions are output.
              image     Do stuff with image data.
           metadata     Extract the ScanImage metadata section from the file.
               help     Print detailed help for the specified sub-command.

Sub-commands
------------

Each sub-command has it's own arguments.  The available sub-commands are listed
below.

.. program:: ScanImageTiffReader

.. option:: help <subcommand>

    Returns more detailed information for each sub-command.

.. option:: descriptions [--frame <index>] <input tiff> [<output file>]

    Extract header data.  :program:`ScanImage` stores (a lot) of information about
    the state of the microscope in each frame.  This command will extract
    all the information for all the frames.  If you use :option:`--frame`
    only that frame's header will be output.

.. option:: metadata <input tiff> [<output file>]

    Extracts the :ref:`metadata` section.  :program:`ScanImage` stores some
    non-frame varying data in the Tiff file, which is what this command
    extracts.

    Historically, this data has had several different formats.  Some are
    human readible text, but some are binary formats that are difficult
    to interpret outside of Matlab.  :program:`ScanImage` provides a
    matlab-specific tool for decoding the binary formats.

.. option:: image bytes <input file>

    Return the size of the raw image data in bytes.

.. option:: image raw <input file> <output file>

    Write raw bytes to an output file for the entire volume.  The pixels are
    written in row-major order (C-style), so for an array with dimensions
    ``[w,h,d]`` the index of a pixel at ``(x,y,z)`` is ``x+y*w+z*w*h``.

.. option:: image shape <input file>

    Print the dimensions and pixel type of the volume in the tiff stack.

.. _ScanImage: http://scanimage.org
