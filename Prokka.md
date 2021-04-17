## Installation

### Bioconda
If you use [Conda](https://conda.io/docs/install/quick.html)
you can use the [Bioconda channel](https://bioconda.github.io/):
```
conda install -c conda-forge -c bioconda -c defaults prokka
```

### Brew
If you are using the [MacOS Brew](http://brew.sh/) 
or [LinuxBrew](http://brew.sh/linuxbrew/) packaging system:
```
brew install brewsci/bio/prokka
```

### Docker
Maintained by https://hub.docker.com/u/staphb
```

docker pull staphb/prokka:latest
docker run staphb/prokka:latest prokka -h
```

### Singularity
```
singularity build prokka.sif docker://staphb/prokka:latest
singularity exec prokka.sif prokka -h
```

### Ubuntu/Debian/Mint
```
sudo apt-get install libdatetime-perl libxml-simple-perl libdigest-md5-perl git default-jre bioperl
sudo cpan Bio::Perl
git clone https://github.com/tseemann/prokka.git $HOME/prokka
$HOME/prokka/bin/prokka --setupdb
```

### Centos/Fedora/RHEL
```
sudo yum install git perl-Time-Piece perl-XML-Simple perl-Digest-MD5 perl-App-cpanminus git java perl-CPAN perl-Module-Build
sudo cpanm Bio::Perl
git clone https://github.com/tseemann/prokka.git $HOME/prokka
$HOME/prokka/bin/prokka --setupdb
```

### MacOS
```
sudo cpan Time::Piece XML::Simple Digest::MD5 Bio::Perl
git clone https://github.com/tseemann/prokka.git $HOME/prokka
$HOME/prokka/bin/prokka --setupdb
```

## Test

* Type `prokka` and it should output its help screen.
* Type `prokka --version` and you should see an output like `prokka 1.x`
* Type `prokka --listdb` and it will show you what databases it has installed to use.
