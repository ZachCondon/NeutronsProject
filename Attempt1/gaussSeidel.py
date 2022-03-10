# Gauss Seidel function


def gaussSeidel(guess,F,S,i,F_size):
    first_sum = 0
    for j in range(0,i):
        first_sum = first_sum + F[i][j]*guess[j]
    second_sum = 0
    for j in range(i+1,F_size):
        second_sum = second_sum + F[i][j]*guess[j]
    return (1/F[i][i])*(S[i]-first_sum-second_sum)
