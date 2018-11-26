#/bin/bash
#not much of an installation. But softlink generation
#must be executed as sudo
currdir=$(pwd)
chmod +x feature_extract.py
sudo ln -s $currdir/feature_extract.py /usr/bin/ST_feature_extract
echo "Installation Complete"

