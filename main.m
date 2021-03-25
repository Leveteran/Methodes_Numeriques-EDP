%
% Resolution numerique du probleme de diffusion avec un 
% schema de Differences Finies, defini sur un intervalle
% ------------------------------------------------------
% -d^2(u)/dx^2 + gamma u    = f in \Omega = ]a, b[
%                      u(a) = g_a
%                      u(b) = g_b
% ------------------------------------------------------
%
% gamma depend de la variable en espace
% -------------------------------------------------


clear all
format long e

disp('Debut du programme')

a = 0 ;
b = 1 ;

% Donnees au bord du domaine 

g_a = 0 ;
g_b = 0 ; 

% Appel de la fonction traitant le probleme

dirichlet1d(a, b, g_a, g_b) ;

disp('Fin du programme') 

