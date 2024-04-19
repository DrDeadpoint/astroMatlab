function newProblemSetup = rewriteProblem(problemSetup,soln)
% newProblemSetup = rewriteProblem(problemSetup,soln)
% used for continuation: plugs in solution values to segments in the problemSetup
% user then replaces values used in continuation scheme
x = soln.x;
solverOptions = problemSetup.options.solverOptions;
fcnOptions = problemSetup.options.fcnOptions;
param = problemSetup.param;
[segs,nodes] = getSegsNodes(x,param,fcnOptions);
param.segments = segs;
param.nodes = nodes;
newProblemSetup = buildProblem(param,solverOptions,fcnOptions);
end