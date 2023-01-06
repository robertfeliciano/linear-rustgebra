use std::fs;

#[derive(Debug)]
pub struct Matrix {
    pub rows: u32,
    pub cols: u32,
    pub data: Vec<Vec<f64>>,
}

impl Matrix {
    pub fn new(rows: u32, cols: u32) -> Matrix {
        let data = vec![vec![0.0; cols as usize]; rows as usize];
        return Matrix { rows, cols, data };
    }

    /** 
     * constructs a matrix from a file
     * example file: 
     * 1.0 0.2 2.0
     * 3.0 4.5 1.2
     * 9.8 3.3 1.4
     */
    pub fn from_file(path: &str) -> Matrix {
        let content = fs::read_to_string(path).unwrap_or_else(|e| panic!("{e}"));
        let mut matrix: Vec<Vec<f64>> = Vec::new();
        for rows in content.lines() {
            let mut row: Vec<f64> = Vec::new();
            let entries: Vec<&str> = rows
                                    .split_whitespace()
                                    .collect();
            for ent in entries {
                row.push(ent.parse::<f64>().unwrap());
            }
            matrix.push(row);
        }

        let r = matrix.len() as u32;
        let c = matrix[0].len() as u32;
        return Matrix { rows: r, cols: c, data: matrix };
    }

    pub fn from_string(input: &str) -> Matrix {
        let mut data: Vec<Vec<f64>> = Vec::new();
        let rows: Vec<&str> = input
                              .split(";")
                              .collect();

        for row in rows {
            let entries: Vec<&str> = row
                                    .split_whitespace()
                                    .collect();
            let mut tmp_row: Vec<f64> = Vec::new();
            for ent in entries {
                tmp_row.push(ent.parse::<f64>().unwrap());
            }
            data.push(tmp_row);
        }
        let n_r = data.len() as u32;
        let n_c = data[0].len() as u32;
        return Matrix { rows: n_r, cols: n_c, data };
    }

    /**
     * returns a deep copy of the matrix
     */
    pub fn copy(&self) -> Matrix {
        let mut matrix: Vec<Vec<f64>> = Vec::new();
        for row in &self.data {
            matrix.push(row.to_vec());
        }
        return Matrix { rows: self.rows, cols: self.cols, data: matrix } ;
    }

    /** 
     * prints a matrix
     */
    pub fn print(&self) {
        for row in &self.data {
            println!("{:?}", row);
        }
        println!();
    }

    /** 
     constructs an identiy matrix
     # Examples
     let mut id = Matrix::new(3,3);

     id.identity();
     */
    pub fn identity(&mut self) {
        if self.rows != self.cols {
            panic!("Not a square matrix.");
        }
        let mut r: u32 = 0;
        while r < self.rows {
            for c in 0..self.cols{
                if r == c{
                    self.data[r as usize][c as usize] = 1.0;
                }
                else {
                    self.data[r as usize][c as usize] = 0.0;
                }
            }
            r += 1;
        }
    }

    /**
     * applies a function over all elements in the matrix
     */
    pub fn apply(&mut self, f: impl Fn(f64) -> f64) {
        for row in 0..self.rows {
            for col in 0..self.cols {
                self.data[row as usize][col as usize] = f(self.data[row as usize][col as usize]);
            }
        }   
    }

    /**
     * adds two matrices and returns the sum
     */
    pub fn add(&self, b: Matrix) -> Matrix {
        if self.rows != b.rows || self.cols != b.cols {
            panic!("Matrices must be of the same size.");
        }
        let mut sum = Matrix::new(self.rows, self.cols);
        for i in 0..self.rows as usize{
            for j in 0..self.cols as usize{
                sum.data[i][j] = self.data[i][j] + b.data[i][j];
            }
        }

        return sum;
    }

    /**
     * subtracts two matrices and returns the difference
     */
    pub fn subtract(&self, b: Matrix) -> Matrix {
        if self.rows != b.rows ||  self.cols != b.cols {
            panic!("Matrices must be of the same size.");
        }
        let mut diff = Matrix::new(self.rows, self.cols);
        for i in 0..self.rows as usize{
            for j in 0..self.cols as usize{
                diff.data[i][j] = self.data[i][j] - b.data[i][j];
            }
        }

        return diff;
    }

    /**
     * computes the dot product of two matrices
     # Examples
     a.dot(b) computes a dot b
     */
    pub fn dot(&self, b: Matrix) -> Matrix {
        if self.rows != b.cols || self.cols != b.rows {
            panic!("Dimensions not matched. M1 is {} by {}, M2 is {} by {}", self.rows, self.cols, b.rows, b.cols);
        }
        let mut dp = Matrix::new(self.rows, b.cols);
        for i in 0..self.rows as usize {
            for j in 0..b.cols as usize {
                let mut sum = 0.0;
                for k in 0..b.rows as usize {
                    sum += self.data[i][k] * b.data[k][j];
                }
                dp.data[i][j] = sum;
            }
        }
        return dp;
    }

    /**
     * computes the Row-Reduced Echelon Form of the matrix
     */
    pub fn rref(&mut self) {
        if self.data[0][0] == 0.0 {
            swap_rows(self, 0);
        }
        let mut lead: usize = 0;
        let rows = self.rows as usize;
        while lead < rows {
            for r in 0..rows{
                let div = self.data[lead][lead];

                let mult = self.data[r][lead] / div;

                if r == lead {
                    //we are at the pivot's row so we need to make it equal to one
                    self.data[lead] = self.data[lead]
                                        .iter()
                                        .map(|entry| entry / div)
                                        .collect();
                }
                else {
                    for c in 0..self.cols as usize {
                        self.data[r][c] -= self.data[lead][c]*mult;
                    }
                }
            }
            lead += 1;
        }
        correct(self);

    }

    /**
     * computes the cofactor of a matrix by expanding upon row "expanded_row"
     */
    pub fn cofactor(&self, expanded_row: u32, j: u32) -> f64 {
        let mut cut: Vec<Vec<f64>> = Vec::new();
        let mut n_r = 0;
        let mut n_c = 0;

        for r in 0..self.rows {
            if r == expanded_row {
                continue;
            }
            n_c = 0;
            let mut v: Vec<f64> = Vec::new();
            for c in 0..self.cols {
                if c == j {
                    continue;
                }
                v.push(self.data[r as usize][c as usize]);
                n_c += 1;
            }
            cut.push(v);
            n_r += 1;
        }

        let minor = Matrix { rows: n_r, cols: n_c, data: cut }.det();
        let base: i32 = -1;
        return minor * f64::from(base.pow(expanded_row + j));
    }

    /**
     *  computes the determinant of a matrix 
     */
    pub fn det(&self) -> f64 {
        if self.rows != self.cols {
            panic!("Determinant requires matrix to be a square matrix. Input matrix was {:?}.", self);
        }
        if self.rows == 2 && self.cols == 2 {
            return self.data[0][0]*self.data[1][1] - self.data[0][1]*self.data[1][0];
        }
        else {
            // matrix is more than 2x2
            // we will just always expand upon row 2 (technically row 1 since the vector is 0-indexed)
            let row: u32 = 1;
            let mut det = 0.0;

            for j in 0..self.data[row as usize].len(){
                det += self.cofactor(row, j as u32) * self.data[row as usize][j];
            }
            return det;
        }
    }

    /**
     * computes the transpose of a matrix
     */
    pub fn transpose(&self) -> Matrix {
        let mut t = Matrix::new(self.cols, self.rows);

        for i in 0..self.rows as usize{
            for j in 0..self.cols as usize{
                t.data[j][i] = self.data[i][j];
            }
        }

        return t;
    }

    pub fn inverse(&self) -> Matrix {
        let mut inv = Matrix::new(self.rows, self.cols);

        for row in 0..self.rows as usize{
            for col in 0..self.cols as usize{
                inv.data[row][col] = self.cofactor(row as u32, col as u32);
            }
        }

        correct(&mut inv);

        inv = inv.transpose();
        let d = self.det();
        if d == 0.0 {
            panic!("Determinant is 0. No inverse.");
        }
        inv.apply(|x| x / d);
        return inv;
    }

    
}

/**
 finds a row with a 1 column and swaps it with the indicated one
 */
fn swap_rows(m: &mut Matrix, row: usize) {
    let mut n_r = 0;
    for r in 0..m.rows as usize{
        if m.data[r][0] > 0.0 {
            n_r = r;
            break;
        }
    }
    let temp: Vec<f64> = m.data[row].clone();
    m.data[row] = m.data[n_r].clone();
    m.data[n_r] = temp;
}

/**
correct negative "zeroes" and floating point approx (1.999999 = 2)
*/
fn correct(m: &mut Matrix) {
    for row in 0..m.rows as usize{
        for col in 0..m.cols as usize{
            let elem = m.data[row][col];
            if elem == -0.0{
                m.data[row][col] = 0.0;
            }
            let floored = elem.floor();
            if elem - floored > 0.9999999{
                m.data[row][col] = elem.round();
            }
            if elem < 0.000001 {
                m.data[row][col] = 0.0;
            }
        }
    }
}