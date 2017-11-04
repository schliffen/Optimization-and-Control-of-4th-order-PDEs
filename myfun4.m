% the objective function
function f = myfun4(Y,M,Bu,ny)

f=.5*(Y(1:ny)-ones(ny,1))'*M*(Y(1:ny)-ones(ny,1))+.05*Y(ny+1,2*ny)'*Bu*Y(ny+1:2*ny);
end
