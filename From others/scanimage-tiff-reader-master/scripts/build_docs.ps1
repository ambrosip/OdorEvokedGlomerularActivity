cd doc
New-Item -ItemType Directory -Force -Path build/html
doxygen
& ./make.bat html
cd ..
Move-Item -Path doc/build/html -Destination public 