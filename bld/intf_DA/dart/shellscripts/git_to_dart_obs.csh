#!/bin/csh -f

cat << EndOfText >! exclude-list.txt
EndOfText

set DEST = $HOME/DART/observations/obs_converters/text_snex/
set SOURCE = $HOME/mpdaf/bld/intf_DA/dart/tools/text_snex/

rsync -Cavzr --exclude-from=exclude-list.txt ${SOURCE} ${DEST}

echo "to finish up ..."
echo "cd ${DEST}"
echo "Go to work directory and do a quickbuild..."
