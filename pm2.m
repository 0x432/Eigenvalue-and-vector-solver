% This function will calcualte a 3x3 symmetric matrixes eigenvalues and
% vectors using SPM2 method and a deflation process.
function [lambda2, x2, B] = pm2(A, lambda1, x0) % Initial start of the PM2.m function takes in the following argumnets A is the matrix, lambda1 is the first eigenvalue, x0 is the initial guess

% The following 4 lines below will check if the matrix is a 3x3 matrix.
% Firstly we need to get the row and column sizes of A we can use the size
% function to get the amount of rows & columns we then store them in a
% variable array. the row will store the amount of rows and the columns will
% store the amount of columns. the if statement line will then check using
% the operators '~=' (not equal to) , '||' (or) if the matrix is the
% correct size the code will carry on if not the code will throw an error to
% end the program if not the rest of the code can execute.
[row, column] = size(A); % Gets the size of the matrix
if row ~= 3 || column ~= 3 % Checks if the rows & columns equal to 3
   error('PM2: The matrix must be a 3x3 matrix.'); % Display an error message if it's not a 3x3 matrix
end % Ends the if statement

% The following 3 lines below will check if the matrix is symmetric. The if
% statement will use the operator '~' which inverts the result of the "is
% equal" function so the statement basically says if not equal. To check if
% the matrix is symmetric we need to plug the original matrix in and then
% the transposed version using the ' operator which transposes the matrix.
% If the transpose of a matrix is not the same as the original matrix,
% then the matrix is not symmetric and we will throw an error to end the
% program. if not the rest of the code can execute.
if ~isequal(A, A') % Checks if the matrix is symmetric
   error('PM2: The matrix must be symmetric.'); % Display an error message if it's not symmetric
end % Ends the if statement

% The following 4 lines below will check if the vector is actually a matrix
% vector. the size function will get the size of the vector and store the
% amount of rows in the rows & the amount of columns in the column section
% of the variable array. Then the if statement will check using the operators
% '~=' (not equal to) , '||' (or) if the vector is the correct size the
% code shall continue to execute if not we shall throw and error and exit
% the execution of the program.
[row, column] = size(x0); % Gets the size of the vector
if row ~= 3 || column ~= 1 % Checks if the vector is the correct size (3x1)
   error('PM2: The vector must be a 3x1') % Displays an error message if it's not a 3x1 (not a vector)
end % Ends the if statement

% The following 2 lines of code deflate matrix A into a new matrix called B
% The purpose of deflation is to remove the first eigenvalue from the
% matrix. Here we use the shifted power method to deflate the matrix. We
% start by creating a varible 'I' which stores an identity matrix the size
% of matrix A. An idtentity matrix is a square matrix with ones on the main
% diagonal and zeros elsewere. To create this we use the built in eye()
% function within matlab. The next line of code starts by creating a
% varible 'B' which is the deflation matrix of A. the following code
% subtracts Lambda1 (first eigenvalue) * The identity matrix from A
% (initial matrix). This gets us our deflated matrix B
I = eye(size(A)); % Creating an Identity marix the size of A
B = A - lambda1 * I; % Deflation process and storeing result in varible B

% Defining Varibles
x2 = x0; % setting x2 as the initial guess
tolerance = 0.00001; % settings the tolerance value 
total_iterations = 100; % Sets the max iterations that can be used
lambda2 = 0; % Settings lambda2 as 0
iteration_counter = 0; % Stores the interations of the while loop

% The following while loop uses SPM2 (Scaled Power Method 2) method to calculate the second eigenvalue and
% vectors of a 3x3 symmetric matrix. The process begins by adding 1 to the
% iteration counter to make sure our loop only runs to the total iteration
% count. Then the processes multiplies the deflated matrix 'B' by x2 which
% is the guess varible and this is stored as 'x'. Next, theo code
% identifies the element in 'x' with the maximum absolute value. This is
% done by taking the absolute values of all elements in 'x', finding the
% maximum among these, and then recording the index of this maximum value.
% The value at this index in 'x' is then assigned to the varible 'm'. This
% step is crucal as 'm' represents the value used for normalizing the
% curernt eigenvector estiamte and updating the eigenvalue estimate in the
% iterative process. the '~' operator basically means ignore this so the
% 'val' is the only thing that is stored within '[]', then x2 is updated by
% multiplying 'x' by m this process is also iterative. Then we create an if
% statment that uses the abs() function within MATLAB which gets the
% absolute value of lambda2 - m then compares with the tolerance and check
% if the difference is smaller than the tolerance if true the break
% function breaks out of the while loop. The next line of code stores m as
% lambda2. 
while iteration_counter < total_iterations % While loop for calculating the second eigenvalue and eigenvector
    iteration_counter = iteration_counter + 1; % Increment the iteration counter

    x = B * x2; % Multiplying the deflated matrix B by the current eigenvector estimate
    [~, val] = max(abs(x)); % Finding the index of the maximum absolute value of 'x', '~' is used to remove the value its self as it is not needed.
    m = x(val); % setting varible 'm' as the element in 'x' with the highest absolute magnitute. 
    x2 = x / m; % Updating the eigenvector estimate
 
    if abs(lambda2 - m) < tolerance % check the difference between the current eigenvalue and the previous eigenvalue to see if its within the tolerance level
        break; % If within the telerance limit break out of the while loop
    end % End of the if statement
    lambda2 = m; % Updating lambda2 to the current value of m
end % End of the while loop

x2 = x2 / norm(x2); % Normalizing the eigenvector to fit with eig(A)

end % End of the PM2 function