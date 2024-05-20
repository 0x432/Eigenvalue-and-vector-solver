% This function will calculate a symmetric 3x3 matrixes eigenvalue and
% eigenvector using Rayleigh quotient and the power method.
function [lambda1, x1] = pm1(A, x0) % Initial start of the PM1 function taking in arguments A & x0 to spit out variables lambda1 (eigenvalue) and x1 (eigenvector)

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
   error('PM1: The matrix must be a 3x3 matrix.'); % Display an error message if it's not a 3x3 matrix
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
   error('PM1: The matrix must be symmetric.'); % Display an error message if it's not symmetric
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
   error('PM1: The vector must be a 3x1') % Displays an error message if it's not a 3x1 (not a vector)
end % Ends the if statement

% Defining Variables
total_iterations = 100; % Total number of iterations the loop will perform
tolerance = 0.00001; % Tolerance variable used to set the threshold of the acceptable error rate
iteration_counter = 0; % Counts and stores how many iterations the loop has gone through
x1 = x0; % x1 initially stores the initial guess

% This while loop uses SMP1 (Scaled Power Method Version 1) To find the
% first eigenvalue and its corresponding eigenvector. The process begins by
% +1 the 'iteration_counter' making sure the loop will only run to the
% 'total_iterations' value. Inside the loop the calculation of the
% x is calculated by multiplying A by x1 which for the first
% iteration will be the initial guess. but will change over iterations.
% The new lambda1 (eigenvalue) is estimated using Rayleigh quotient, which
% is a method that multiplies the inverse of the current eigenvector
% estimate 'x1' by 'x' then divides by the inverse of 'x1'
% * 'x1' uses our current eigenvalue estimate. After this 'x1' is
% normalized by dividing 'x' by the normalized version of 'x' using matlabs
% in built norm() function. This makes sure the magnitude of 'x1' is constant.
% Then and if statement is created to check if the eigenvalue & vector are close
% enough to the actual values within the tolerance of 0.00001.
% If so a break; is used to exit out of the while loop. if the if statement
% is false the loop will continue to iterate to either 100 or until it has
% an eigenvalue / vector estimate that fits the tolerance.
while iteration_counter < total_iterations % While loop for calculating the first eigenvalue and eigenvector
   iteration_counter = iteration_counter + 1; % Increment the iteration counter
  
   x = A * x1; % Multiply matrix A by the current eigenvector estimate x1
   lambda1 = (x1' * x) / (x1' * x1); % Calculate new eigenvalue estimate using Rayleigh quotient
   x1 = x / norm(x); % Normalize the matrix product to get the updated eigenvector estimate x1
   if norm(A * x1 - lambda1 * x1) < tolerance % Check if the current eigenvalue and eigenvector estimate is close enough to the tolerance also called convergence
       break; % Exit the loop if the estimates are within the tolerance level
   end % Ends the if statement
end % Ends the while loop
end % End of of the pm1.m function