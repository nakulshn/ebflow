#!/bin/bash

# generate_figure2_jobs.sh



# Only for block02corr0.9, nbyp1.0, all 4 priors, 10 seeds



for prior in gaussian skew cauchy bimodal; do

  if [ "$prior" = "gaussian" ]; then

    noiser=0.5

    lambda=0.003

  else

    noiser=0.8

    lambda=0.001

  fi

  

  for seed in {1..10}; do

    # EBflow (not preconditioned)

    echo "Rscript ../run_EBflow.R --P $prior --noiser $noiser --lambda $lambda --D block02corr0.9 --dimr 1.0 --eta.phi decay --eta.w decay --seed $seed --respath ../Figure2/results/ --predict"

    

    # EBflow (preconditioned)

    echo "Rscript ../run_EBflow.R --P $prior --noiser $noiser --lambda $lambda --D block02corr0.9 --dimr 1.0 --eta.phi decay --eta.w decay --seed $seed --precondition --respath ../Figure2/results/ --predict"

    

    # Gibbs

    echo "Rscript ../run_Gibbs.R --P $prior --noiser $noiser --lambda $lambda --D block02corr0.9 --dimr 1.0 --T 100 --seed $seed --respath ../Figure2/results/ --predict"

  done

  

  # VI (no seed needed, only once per prior)

  echo "Rscript ../run_VI.R --P $prior --noiser $noiser --lambda $lambda --D block02corr0.9 --dimr 1.0 --respath ../Figure2/results/ --predict"

done
