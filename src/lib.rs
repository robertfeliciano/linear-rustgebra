use std::fs;

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
        let mut c: u32 = 0;
        let mut r: u32 = 0;
        for rows in content.lines() {
            r += 1;
            let mut row: Vec<f64> = Vec::new();
            let entries: Vec<&str> = rows
                                    .split_whitespace()
                                    .collect();
            for ent in entries {
                if r == 1 {
                    c += 1;
                }
                row.push(ent.parse::<f64>().unwrap());
            }
            matrix.push(row);
        }
        return Matrix { rows: r, cols: c, data: matrix };
    }

    pub fn from_string(input: &str) -> Matrix {
        let mut data: Vec<Vec<f64>> = Vec::new();
        let rows: Vec<&str> = input
                              .split(";")
                              .collect();
        let n_r: u32 = rows.capacity() as u32;
        let mut n_c: u32 = 0;

        for row in rows {
            let entries: Vec<&str> = row
                                    .split_whitespace()
                                    .collect();
            let mut tmp_row: Vec<f64> = Vec::new();
            for ent in entries {
                tmp_row.push(ent.parse::<f64>().unwrap());
            }
            n_c = tmp_row.capacity() as u32;
            data.push(tmp_row);
        }

        return Matrix { rows: n_r - 1, cols: n_c - 1, data };
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
        println!("[");
        for row in 0..self.rows {
            print!("   [");
            for col in 0..self.cols {
                if col != self.cols - 1{
                    print!("{} ", self.data[row as usize][col as usize]);
                }
                else {
                    print!("{}", self.data[row as usize][col as usize]);
                }
            }
            println!("]");
        }
        println!("]");
    }

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
    pub fn map(&mut self, f: impl Fn(f64) -> f64) {
        for row in 0..self.rows {
            for col in 0..self.cols {
                self.data[row as usize][col as usize] = f(self.data[row as usize][col as usize]);
            }
        }   
    }

    /**
     * computes the Row-Reduced Echelon Form of the matrix
     */
    pub fn rref(&mut self) {

    }

    /**
     * computes the cofactor of a matrix by expanding upon row "expanded_row"
     */
    pub fn cofactor(&self, expanded_row: u32, j: u32) -> f64 {
        // want to elimate row i and col j from matrix m
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
            // self.print();
            panic!("Determinant requires matrix to be a square matrix. This matrix was {} by {}", self.rows, self.cols);
        }
        if self.rows == 2 && self.cols == 2 {
            return self.data[0][0]*self.data[1][1] - self.data[0][1]*self.data[1][0];
        }
        else {
            // matrix is more than 2x2
            // we will just always expand upon row 2 (technically row 1 since the vector is 0-indexed)
            let row: u32 = 1;
            let mut det = 0.0;

            for j in 0..self.data[row as usize].capacity() - 1 {
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

        // need to fix negative zeroes. idek how that happens tbh
        for row in 0..inv.rows as usize{
            for col in 0..inv.cols as usize{
                if inv.data[row][col] == -0.0{
                    inv.data[row][col] = 0.0;
                }
            }
        }

        inv = inv.transpose();
        let d = self.det();
        if d == 0.0 {
            panic!("Determinant is 0. No inverse.");
        }
        inv.map(&|x| x / d);
        return inv;
    }

    
}
