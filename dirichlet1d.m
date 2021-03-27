function [err, soldir] = dirichlet1d(a, b, g_a, g_b) 

  disp('Choisir un pas de maillage, par exemple h = 0.1 ou 0.01 ou ... ou 0.001 ')

  h = input('Effectuez votre choix, par exemple 0.1 ou 0.01 ou ... ou 0.001  '); 
  disp('Vous avez choisi h = ')
  disp(h)

  % Maillage du domaine 
  % ------------------------------------------------------------

  % Determination du nombre de points de discretisation du domaine

  N = round((b-a)/h);

  h_evaluated = (b-a)/N ;

  % Pas de maillage considere 

  h = h_evaluated ;

  disp('Pas de maillage considere : h = ')
  disp(h)

  % Nombre de grilles 

  N_grilles = N ;

  disp('Nombre de grilles : N_grilles = ')
  disp(N_grilles)

  % Nombre de noeuds internes : Nx_internal
  % Nx_internal = Nombre total de points de discretisation - 2 
  % Nx_internal = (N + 1) - 2 = N - 1

  Nx_internal = N - 1 ;

  disp('Nombre de noeuds internes : Nx_internal = ')
  disp(Nx_internal)

  % Determination des noeuds internes 
  % xdiscrete(i), 1 <= i <= Nx_internal
  % Noter en outre : x_{0} = a, x_{last} = b 

  x0 = a ;

  xdiscrete(1) = x0 + h ;
  for i = 2:Nx_internal 
      xdiscrete(i) = xdiscrete(i-1) + h ;
  end

  xlast = b ; 

  % Determination des coefficients relatifs a chaque noeud interne
  % --------------------------------------------------------------

  % Coefficients provenant de f
  % fvalues_i : valeur de f au ieme noeud 

  for i = 1:Nx_internal 
      fvalues(i) = feval('f', xdiscrete(i)) ;
  end 


  % Saisie de la matrice du systeme 
  % --------------------------------

  A(1:Nx_internal,1:Nx_internal) = 0.0;

  for i = 1:Nx_internal 

      gamma_i  = feval('gammaExpr', xdiscrete(i)) ;      

      A(i,i) = 2 + (gamma_i)*h^2 ;
  end 

  for i = 1:(Nx_internal-1)

      A(i+1,i) = -1 ;

  end 

  for i = 1:(Nx_internal-1) 

      A(i,i+1) = -1 ;

  end

  % Repartition des coefficients non nuls de la matrice du systeme

  disp('Visualisez la repartition des coefficients non nuls de la matrice')
  disp('Et saisir << return >> pour la suite des calculs')
  spy(A)

  keyboard

  % Composantes du second membre du systeme 

  B(1) = (h^2)*fvalues(1) + g_a ; 

  for i = 2:(Nx_internal-1) 

      B(i) = (h^2)*fvalues(i) ; 

  end 

  B(Nx_internal) = (h^2)*fvalues(Nx_internal) + g_b ;

% Calcul de la solution de A X = B :
  % ---------------------------------

  % X : approximation par Differences Finies de (u(xdiscrete(i)))_{1 <= i <= Nx_internal} 

  X = A\B' ;

  % Consideration de la solution discrete sur tout le domaine 

  soldir = [g_a ; X ; g_b] ;

  % Consideration de la solution exacte dans la base discrete 

  for i = 1:Nx_internal
      u(i) = feval('uex', xdiscrete(i)) ;
  end

  u = [g_a, u, g_b] ;

  erreur = norm(u' - soldir, inf)/norm(u', inf) ;
  err    = 100 * erreur ;

  disp('Erreur commise en % avec la methode de discretisation est de :')
  disp(err)

  % Representation graphique de la solution

  x = [x0, xdiscrete, xlast];

  plot(x, soldir, 'bo')
  hold on
  fplot('uex', [a b],'r-')
  legend('Solution discrete', 'Solution exacte', -1)
  xlabel('Abscisses')
  ylabel('Valeurs prises par la solution discrete et la solution exacte')
  title('Representation de la solution discrete et de la solution exacte')

  %Impression en postscript
  print dirichlet1d_10.ps

  hold off

end 

