function[sigma_tr, p_tr, s_tr, rho_tr, n_tr]=compute_trials(eps_tot, eps_p_prev, mu, lambda)

eps_el = eps_tot-eps_p_prev;
trace_el = eps_el(1) + eps_el(4) + eps_el(5);

Id = [1,0,0,1,1]';
VOL = Id*Id';
DEV = diag([1,1,1,1,1]') - VOL/3;

s_tr = DEV*eps_el;

sigma_xx=2*mu*eps_el(1)+lambda*trace_el;
sigma_xy=2*mu*eps_el(2);
sigma_yy=2*mu*eps_el(4)+lambda*trace_el;
sigma_zz=2*mu*eps_el(5)+lambda*trace_el;

sigma_tr=[sigma_xx, sigma_xy, sigma_xy, sigma_yy, sigma_zz]';

p_tr=sigma_xx+sigma_yy+sigma_zz;
            
norm_e_dev = sqrt(max(0,eps_el'*s_tr));

rho_tr=2*mu*norm_e_dev;

 if(abs(rho_tr)<1e-16)              
     disp('Error: rho_tr==0')
 end

n_tr=s_tr/norm_e_dev;