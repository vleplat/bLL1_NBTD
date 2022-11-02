% Ce code vient générer un jeu de données admettant un modèle LL1 exact. Il
% est construit comme spécifié dans
% https://hal.archives-ouvertes.fr/hal-03158076v2/document, Section 6.2.2

clear all
close all
clc

%% Dimensions

R = 4; L = [3 3 6 6]; 
I = 120; J = 120; K = 173;

%% Load spectra

load('end4.mat') %Spectres tirés de Jasper Ridge
C1 = M(1:K,1); C2 = M(1:K,2); C3 = M(1:K,3); C4 = M(1:K,4); 

C = [C1 C2 C3 C4];
clear A cood M

%% Generate abundance maps

A1 = zeros(I,J); A1(1:40,1:40) = 1; A1(41:80,41:80) = 1; A1(81:I,81:J) = 1; 

A2 = zeros(I,J); A2(1:40,41:80) = 1; A2(41:80,81:J) = 1; A2(81:I,1:40) = 1;

A3 = zeros(I,J); A3(41:60,1:20) = 1; A3(61:80,21:40) = 1; A3(101:I,41:60) = 1;
                 A3(81:100,61:80) = 1; A3(21:40,81:100) = 1; A3(1:20,101:J) = 1;
                 
A4 = ones(I,J) - (A1+A2+A3);

S = [A1(:) A2(:) A3(:) A4(:)];

%% Generate SRI 

% for r=1:R
%     C(:,r) = C(:,r)/norm(C(:,r));
%     S(:,r) = S(:,r)/norm(S(:,r));
% end

Z = reshape(S*C',[I J K]);


% figure(1)
% for r=1:R
%     subplot(1,3,r)
%     plot(C(:,r));
% end
% figure(2)
% for r=1:R
%     subplot(1,3,r)
%     imagesc(reshape(S(:,r),[I,J]));
% end
