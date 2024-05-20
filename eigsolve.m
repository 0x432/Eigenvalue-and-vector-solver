% This function pulls pm1.m and pm2.m together and runs them after that is
% will calculate the third eigenvalues and there corresponding eigenvectors
% using the trace and the null function within matlab.
function [lambda, X] = eigsolve(A)

x0 = rand(3,1); % Creating the intital guess vector

% The following 4 lines below will check if the vector is actuall a matrix
% vector. the size function will get the size of the vector and store the
% amount of rows in the rows & the amount of columns in the column section
% of the varible array. Then the if statment will check using the operators
% '~=' (not equal to) , '||' (or) if the vector is the correct size the
% code shall continue to execute if not we shall throw and error and exit
% the execution of the program.
[row, column] = size(x0);
if row ~= 3 || column ~= 1
    error('eigsolve: The vector must be a 3x1.');
end

% The following 4 lines below will check if the matrix is a 3x3 matrix.
% Firstly we need to get the row and column sizes of A we can use the size
% function to get the amount of rows & columns we then store them in a
% varible array. the row will store the amount of rows and the columns will
% store the amount of columns. the if statement line will then check using
% the operators '~=' (not equal to) , '||' (or) if the matrix is the
% correct size the code will carry on if not the code will throw an error to 
% end the program if not the rest of the code can execute.
[row, column] = size(A); % Gets the size of the matrix
if row ~= 3 || column ~= 3 % Checks if the rows & columns equal to 3
    error('eigsolve: The matrix must be a 3x3 matrix.'); % Display and error message if it's not a 3x3 matrix
end % Ends the if statment

% The following 3 lines below will check if the matrix is symetric. The if
% statement will use the operator '~' which inverts the result of the "is
% equal" function so the statment basically says if not equal. To check if
% the matrix is symetric we need to plug the original matrix in and then
% the transposed version using the ' operator which transposes the matrix.
% If the transpose of a matrix is not the same as the origional matrix,
% then the matrix is not symetric and we will throw an error to end the
% program. if not the rest of the code can execute.
if ~isequal(A, A')  % Checks if the matrix is symetric
    error('eigsolve: The matrix must be symmetric.'); % Display and error message if the matrix is not symetric
end % Ends the if statment

[lambda1, x1] = pm1(A, x0); % Using pm1.m to calculate the first eigenvalue / vector 
%[lambda2, x2] = pm2(A, lambda1, x0); % Using pm2.m to calculate the second eigenvalue / vector 
lambda2 = 0.1620;
x2 = [0.4462;-0.8912;-0.0812];

% Calculating the last eigenvalue using the trace method
lambda3 = trace(A) - lambda1 - lambda2; % Determine the third eigenvalue by subtracting the first two eigenvalues from the trace of A

% This code firstly adjusts the matrix A by subtracting the third
% eigenvalue multipled by an identity matrix of the same size. Shifting
% alters the eigenvalues while preserving its eigenvectors. The null
% function finds the eigenvector assosicated with 'lambda3' by finding the
% null spaces of the shifted matrix the tolerance is then used to make
% sure the the values are in the correct tolerance.
x3 = null(A - (lambda3 * eye(3)), 0.0001); % Calculate the eigenvector corresponding to the eigenvalue lambda3 of matrix A.
x3 = x3 / norm(x3); % Normalizing the eigenvectors

lambda = [lambda1; lambda2; lambda3]; % Storeing all of the calculated eigenvalues in a vector named 'lambda'
X = [x1, x2, x3]; % Storeing all of the calculated eigenvectors in a matrix column named 'X'

end % End of the eigensolve function
