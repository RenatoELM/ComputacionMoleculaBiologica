#include <iostream>
using namespace std;

int n;
string A;
int dp[100][100];

bool emp(char X, char Y) {
    return ((X == 'A' && Y == 'U') || (X == 'U' && Y == 'A') || (X == 'C' && Y == 'G') || (X == 'G' && Y == 'C'));
}

int n_j(int i, int j) {
    if (dp[i][j] != -1) {
        return dp[i][j];
    } 

    if (i >= j - 1) {
        dp[i][j] = 0;
        return 0;
    }

    int ret = 0;
    ret = max(ret, n_j(i + 1, j));
    ret = max(ret, n_j(i, j - 1));
    if (emp(A[i], A[j])) {
        ret = max(ret, n_j(i + 1, j - 1) + 1);
    }
    for (int k = i + 1; k < j; k++) {
        ret = max(ret, n_j(i, k) + n_j(k + 1, j));
    }

    dp[i][j] = ret;
    return ret;
}

inline int n_j() {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dp[i][j] = -1;
        }
    }

    return n_j(0, n - 1);
}

string cad(int i, int j) {
    if (i == j) {
        return A.substr(i, 1);
    }
    if (i > j) {
        return "";
    }

    if (n_j(i, j) == n_j(i + 1, j)) {
        string left = A.substr(i, 1);
        return left + cad(i + 1, j);
    }

    if (n_j(i, j) == n_j(i, j - 1)) {
        string right = A.substr(j, 1);
        return cad(i, j - 1) + right;
    }

    if (emp(A[i], A[j]) && n_j(i, j) == n_j(i + 1, j - 1) + 1) {
        string left = A.substr(i, 1);
        string right = A.substr(j, 1);
        return "(" + left + cad(i + 1, j - 1) + right + ")";
    }

    for (int k = i + 1; k < j; k++) {
        if (n_j(i, j) == n_j(i, k) + n_j(k + 1, j)) {
            return cad(i, k) + cad(k + 1, j);
        }
    }
}

int main() {
    A = "GGAAAUCC";
    //A = "ACUCGAUUCCGAG";
    //A = "ACCUAGGAAACUUAGGAUCCAUU";
    n = size(A);

    printf("%d\n", n_j());
    cout << cad(0, n - 1) << endl;

    return 0;
}
