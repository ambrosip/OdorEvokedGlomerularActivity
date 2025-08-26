$root=Split-Path -parent $MyInvocation.MyCommand.Path
$src= join-path $root ../doc
write-host $src

$watcher=new-object System.IO.FileSystemWatcher
$watcher.Path = $src
$watcher.Filter = '*.*'
$watcher.IncludeSubDirectories=$true
$watcher.EnableRaisingEvents=$true
$watcher.NotifyFilter="LastAccess,LastWrite,FileName,DirectoryName"
write-output $watcher

while(1) {
    $watcher.WaitForChanged("All")
    . (join-path $root "build_docs.ps1")
    write-host "---------------------------"
}
