optimizer = WaiterMotionOptimizer(6,2);
q0 = [1.56; -0.57; -0.20; 0; -0.36; 0];
qe = [-1.56; 0.57; 0.20; 0; -0.36; 0];
[X_opt, U_opt, T_opt, stats] = optimizer.solve_optimization(q0, qe);
disp(stats.return_status);
optimizer.plot_results(X_opt, U_opt, T_opt);
