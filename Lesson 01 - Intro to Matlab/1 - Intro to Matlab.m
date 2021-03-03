%% Comments and shortcuts

% if you use a percentage sign, then you can write
% whatever you want ... % means a comment

% Matlab as a calculator
3+4

% use ctrl + R to comment out several lines (you have to select the lines)
5*2
% % 54/3
% % 234^2
144^0.5 

% and use ctrl+T to remove the comment mode
sqrt(144)

%{
54/3
234^2
144^0.5 
%}

%% Defining variables

% use the = sign to assign variables
clc

% use ; to suppress the output from the command window
car = 3;
nameyouwant = 10;
% Matlab is case sensitive
Car = 5

Car*nameyouwant


%% Tic & Toc

tic
var1 = 13;
var2 = 14;
var3 = 234;
toc

% use clear to delete all variables from the workspace

% or use clear variablename to delete one

clear car

clear var*

clear *1

%% Data types

time = 10;

% defining an array of characters

varA = 'souumapalavra';
varB = 'carro';

% defining a string
varC = "agorasimstring"

% from character vector to string
 varD = string(varB)
% string indexing

% from string to character vector
varE = char(varC)

length(varA)

varA(13)


% from position 3 to 8, one by one
varA(3:9)

% from 3 to 9 taking every other
varA(3:2:9)


% Concatenation

varF = [varA varB]

% % stringarray = []

% knowing the respective numbers in the ASCII table

varA+0

%% Arrays

v = [3 19 24]

size(v)

vtranspose = v'

M = [2, 3, 2, 234; 12, 3214, 12 , 1]

% random number from 0 to 1
M2 = rand(5,5)

%% Array indexing

M2

% accessing data via its indices
M2(3,2)

% getting an entire row
M2(3,:)

% a column
M2(:,2)

M2(:,[2 4])

M2(:,[4 2])

% you can use vector for indexing

v = [1 3]
M2(:,v)

v3 = 3:4:23

% use end to access the last entry
v3(end)

M2([3,4],3:end)

%% You can use indexing to change values

M3 = eye(5)

M3(1,3) = 43

M3(end,2:4) = 100

% deleting columns/rows

M4 = eye(6)

M4(:,[2, 4]) = []

M4(2,:) = []

%% Defining a vector through the function linspace

v3 = linspace(1,5,200)

%% Matrix operations

M1 = [1 3; 4 5]
M2 = [0 14;-5 pi]

M1+M2
M1-M2

% matricial multiplication
M1*M2

% point by point multiplication
M1.*M2

% operation between a matrix and a scalar
3*M1
M1/10


%% Logical operators

answer = 3>300

answer = 300==300 % equality
answer = 3~=300 % different




M5 = [-1 234 42; 10 50 243; -20 -236 90]

I = M5>0

% you can use logical operators for indexing

% M5(I) = 1000

M5(M5>0) = 1000


%% AND and OR operators

clc

condition = 3>0 & 5>30 
condition = 3>0 | 5>30 

v = randn(1,10)

%% Function Find

v > 0

I = find(v>0)
v(I)


%% IF, ELSE, ELSEIF
clc
word = 'hello world';
if 300>20
    disp(word)
elseif 300>30
    sqrt(100)
else
    3+4
end


%% For and While loops

clc
% vector = 1:2:8
vector = [0.2342 0.763 1.234]
for j= vector
%     disp(['j equals ' int2str(j)])
 j
% 3+4
end
    
%% Loops inside loops
clc
clear vv
for j=1:5
    vv(j) = 2*j;
    for i=1:5
   MM(i,j) = i*j;   
    end
end

%% While

count = 0

while count < 10
    disp('lower than 10')
    count = count+1
end

%% Cells

C0{4,4} = 4

C{1,1} = 1;
C{1,2} = 'hello world';
C{2,1} = C0;
C{2,2} = randn(3,5);


%% Structs

S = dir

%%

exp(1).ratname = 'mickeymouse';
exp(1).data = 'May';
exp(1).EEG = randn(1,1000);

exp(2).ratname = 'minie';
exp(2).data = 'Jun';
exp(2).EEG = [1 2 3 4 523 5213 1];

%% Saving data

% save everything
save('lesson1.mat')

% save some variables
save('lesson1b.mat','C0','M1','S')

% appending variables to an existing file

filename = 'lesson1b.mat'
save(filename,'var*','Car','-append')

%% Loading data

clear

% load everything
% load('lesson1.mat')

load('lesson1.mat','var*')




    
    






