# Install Circos

1. Download *find link*
2. Extract
3. Install required Perl modules (install libgd for the GD module in perl)

```{bash}
sudo pacman -Syu gd
chmod +x bin/circod
sudo cpan $(bin/circos -modules | grep "^missing" | rev |cut -d' ' -f1 | rev | sed -z 's/\n/ /g')
bin/circos -modules

wget http://www.cpan.org/authors/id/L/LD/LDS/GD-2.49.tar.gz
tar -xzvf GD-2.49.tar.gz
cd GD-2.49
perl Makefile.PL
make
sudo make install

sudo cpan GD
sudo cpan GD::Polyline
```
