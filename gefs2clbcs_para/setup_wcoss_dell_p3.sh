module purge

module use /usrx/local/nceplibs/dev/hpc-stack/libs/hpc-stack/modulefiles/stack
module load hpc/1.1.0
module load hpc-ips/18.0.1.163
module load hpc-impi/18.0.1
module load bacio/2.4.1
module load nemsio/2.5.2
module load w3emc/2.7.3
module load w3nco/2.4.1

module use /gpfs/dell2/emc/modeling/noscrub/emc.nemspara/soft/modulefiles
module load hdf5_parallel/1.10.6
module load netcdf_parallel/4.7.4

#echo $NETCDF

export NETCDF_LDFLAGS="-L$NETCDF/lib -lnetcdff -lnetcdf"
