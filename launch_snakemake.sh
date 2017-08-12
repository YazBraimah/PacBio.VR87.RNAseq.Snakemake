if [ "$#" -lt 1 ]; then
    echo "script.sh <out_dir (with config.yaml)> [additional snakemake arguments]*"
    exit
fi

outdir=$1

mkdir -p $outdir/cluster_logs

snakemake -j 16 --local-cores 4 -w 90 --cluster-config cluster.json --cluster "qsub -S /bin/bash -q {cluster.q} -j y -o {cluster.output} -N {cluster.jname} -m be -pe bscb {cluster.cores} -binding linear:{cluster.bl} -l h_vmem={cluster.memory} -l h_rt={cluster.time}" --directory "$@"
