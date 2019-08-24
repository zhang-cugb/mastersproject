# Setup computer for running PorePy
1. Download the latest version of VirtualBox.
2. Download the latest version of Ubuntu (https://ubuntu.com/download/desktop)

### Setup a new virtual machine
Specifications: 
1. Ubuntu 64-bit (if 64-bit not found, check out BIOS settings for virtualization and VT-x)
2. Base memory: 4096MB, Storage >=25GB
3. Minimal ubuntu installation

After ubuntu installation, navigate to "Devices" at top of VM window, click "Insert Guest Additions CD image....". This will fix issues with screen resolution

### Prepare PorePy install
Update system and install basic features:
```sudo apt update && yes|sudo apt upgrade && sudo apt install -y wget git```

Setup a dedicated folder for all porepy related files
```
cd ~/ 
mkdir porepy
cd porepy
```

Install `gmsh` for meshing.
```
wget http://gmsh.info/bin/Linux/gmsh-4.3.0-Linux64.tgz
sudo tar xvf gmsh-4.3.0-Linux64.tgz
sudo rm gmsh-4.3.0-Linux64.tgz
echo "config={\"gmsh_path\":\"~/porepy/gmsh-4.3.0-Linux64/bin/gmsh\"}" > porepy_config.py
```

Install point_in_polyhedron script.
```
cd ~/porepy
wget https://raw.githubusercontent.com/keileg/polyhedron/master/polyhedron.py
mv polyhedron.py robust_point_in_polyhedron.py
```

Install miniconda3
```
cd /tmp
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O conda.sh
chmod +x conda.sh
bash conda.sh -bp ~/.conda
export PYTHONPATH=~/porepy
```
Note: select "yes" if asked to initialize conda, then restart shell.

Update conda
```
conda update conda
```
### Install PorePy:
Make sure you are in (base) conda env.
```
cd ~/porepy
git clone https://github.com/pmgbergen/porepy.git pp
pip install -r requirements-dev.txt
cd pp
pip install .
```

Test installation:
```
cd ~/porepy/pp
python setup.py test
```

You should get something like 934 passed. I got 80% total cover.
The start of my test session looked like this:
```
====================================== test session starts ======================================
platform linux -- Python 3.7.3, pytest-4.5.0, py-1.8.0, pluggy-0.12.0
rootdir: /home/haakon/porepy, inifile: setup.cfg
plugins: cov-2.6.1
collected 934 items
```

