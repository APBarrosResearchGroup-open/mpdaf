#!/bin/csh -f

cat << EndOfText >! exclude-list.txt
EndOfText

set DEST = $HOME/DART/models/mshm/
set SOURCE = $HOME/mpdaf/bld/intf_DA/dart/mshm/

rsync -Cavzr --exclude-from=exclude-list.txt ${SOURCE} ${DEST}

echo "to finish up ..."
echo "cd ${DEST}"
echo "Go to work directory and do a quickbuild..."
