% Max tolerated connection divergence
mu=0.2;
%par_file_name=[tempname('/tmp') '-%s.raw']
par_file_name='/tmp/par-%s.raw';
rand_par_file_sin_noise(25,
                        20,
                        par_file_name,
                        1,
                        0.01,
                        0.1,
                        0.1,
                        1,
                        10);
par_file_name=sprintf(par_file_name,'par');
% Connection analysis Mcaulay & Quatieri
%mq_cxns_path=[tempname('/tmp') '-%s.raw']
mq_cxns_path='/tmp/mq-%s.raw';
par_cxn_mq(par_file_name,mq_cxns_path,mu);

% Connection analysis LP with simple divergence cost
%lp_cxns_path=[tempname('/tmp') '-%s.raw']
lp_cxns_path='/tmp/lp-%s.raw';
[costs_lp,ernums_lp]=par_cxn_lp_b(par_file_name,
             lp_cxns_path,
             1,
             0,
             2,
             0.95,
             mu);

% Connection analysis LP with cost based on prediction error
%lp_cxns_b_path=[tempname('/tmp') '-%s.raw']
lp_cxns_b_path='/tmp/lp_b-%s.raw';
[costs_lp_b,errnums_lp_b]=par_cxn_lp_b(par_file_name,
             lp_cxns_b_path,
             0,
             1,
             2,
             0.95,
             mu);

% Connection analysis LP with cost based on prediction error
% and simple distance criterion
lp_cxns_ab_path='/tmp/lp_ab-%s.raw';
[costs_lp_ab,errnums_lp_ab]=par_cxn_lp_b(par_file_name,
             lp_cxns_ab_path,
             0.5,
             0.5,
             2,
             0.95,
             mu);

cxns_mq=plot_par_cxn_files(par_file_name,mq_cxns_path,0,1);
title('Mcaulay & Quatieri');
cxns_lp=plot_par_cxn_files(par_file_name,lp_cxns_path,0,2);
title('LP simple cost');
cxns_lp_b=plot_par_cxn_files(par_file_name,lp_cxns_b_path,0,3);
title('LP prediction error cost');
cxns_lp_ab=plot_par_cxn_files(par_file_name,lp_cxns_ab_path,0,4);
title('LP prediction error and simple cost');
