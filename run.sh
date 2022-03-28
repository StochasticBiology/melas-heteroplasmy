commandstr=$1

if [ $# -eq 0 ]
then
    echo "Running demo by default: use --full for full analysis."
    commandstr="--demo"
fi

if [[ $commandstr == *--demo* ]]; then
  cd src/simulate
  gcc -o3 met_hast.c -lm -o met_hast.ce -g
  ./met_hast.ce 1 0.1 42
  cd ../plots
  python3 plot_posterior_and_MAP.py 1
  python3 plot_mcmc_extension.py 1
  python3 fancy_plot_mcmc_extension.py 1
fi

if [[ $commandstr == *--full* ]]; then
  cd src/simulate
  gcc -o3 met_hast.c -lm -o met_hast.ce -g
  ./met_hast.ce 0 0.1 42
  cd ../plots
  python3 plot_posterior_and_MAP.py 0
  python3 plot_mcmc_extension.py 0
  python3 fancy_plot_mcmc_extension.py 0
fi
