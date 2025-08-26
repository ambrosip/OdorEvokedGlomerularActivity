function Get-Batchfile ($file) {
    $cmd = "`"$file`" amd64 & set"
    cmd /c "$cmd" | Foreach-Object {
        $p, $v = $_.split('=')
        Set-Item -path env:$p -value $v
    }
}
  
function VsVars32([int]$version = 10)
{
    $versionkey = $version.ToString("00.0")

    $batfilename = "vcvarsall.bat" 

    $key = "HKCU:SOFTWARE\Microsoft\VisualStudio\" + $versionkey + "_Config"

    write-host $key

    $VsKey = get-ItemProperty $key
    $VsInstallPath = [System.IO.Path]::GetDirectoryName($VsKey.ShellFolder)
    write-host $VsInstallPath
    $VsToolsDir = [System.IO.Path]::Combine($VsInstallPath, "VC")
    $BatchPath = join-path $VsToolsDir $batfilename
    
    Get-Batchfile( $BatchPath )

    write-host (Get-Command cl).path

    [System.Console]::Title = "PS MSVS" + $versionkey
}

cd build
VsVars32(14)
msbuild ScanImageTiffReader.sln /property:Configuration=Release
cd ..

