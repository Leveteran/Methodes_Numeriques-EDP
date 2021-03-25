function [fd] = f(x)

  % The function f
  % (depending below on gamma)

  % Consideration of the function gamma

  gammafunc = @gammaExpr;


  % Case where u = x(1-x) 
  % ---------------------
  %
  % f = 2 + gamma*(x(1-x))  

%  fd = 2 + gammafunc(x).*(x.*(1-x)); 


  % Case where u = sin(pi*x)
  % ------------------------
  %
  % f = pi^2*sin(pi*x) + gamma*sin(pi*x) 

  fd = pi^2*sin(pi*x) + gammafunc(x).*sin(pi*x); 

end

