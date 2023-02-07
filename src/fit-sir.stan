functions {
  vector sir(real t, vector y, real beta, real alpha) {
    vector[2] dydt;
    dydt[1] = -beta * y[1] * y[2];
    dydt[2] = beta * y[1] * y[2] - alpha * y[2]; // infected update
    return dydt;
  }
}

data {
    int max_t; // duration of simulation, in days
    int max_obs_t;
    real ts[max_obs_t+1]; // time points to evaluate SIR
    int<lower=0> y[max_t+1]; // cases for each day
    real I0;
}

parameters {
    real<lower=0.3, upper=1.5> beta;
    real<lower=0.05, upper=0.85> alpha; 
    real<lower=0.1, upper=1-I0> S0;
}

transformed parameters {
    real<lower=0, upper=1> inf_curve[max_obs_t+1];
    {
        vector[2] init;
        init[1] = S0;
        init[2] = I0;
        inf_curve = ode_rk45(sir, init, -0.0001, ts, beta, alpha)[:,2];
    }
}

model {
    for (t in 1:(max_obs_t+1))
      y[t] ~ poisson(1000 * inf_curve[t]);
}
