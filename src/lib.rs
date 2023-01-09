use std::fs;

#[derive(Debug)]
pub struct Matrix {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<Vec<f64>>,
}

impl Matrix {
    pub fn new(rows: usize, cols: usize) -> Matrix {
        let data = vec![vec![0.0; cols]; rows];
        return Matrix { rows, cols, data };
    }

    pub fn from_file(path: &str) -> Matrix {
        let content = fs::read_to_string(path).unwrap_or_else(|e| panic!("{e}"));
        let mut matrix: Vec<Vec<f64>> = Vec::new();
        for rows in content.lines() {
            let mut row: Vec<f64> = Vec::new();
            let entries: Vec<&str> = rows
                                        .split_whitespace()
                                        .collect();
            
            entries.iter().for_each(|ent| row.push(ent.parse::<f64>().unwrap()));

            matrix.push(row);
        }
        let r = matrix.len();
        let c = matrix[0].len();
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

            entries.iter().for_each(|ent| tmp_row.push(ent.parse::<f64>().unwrap()));

            data.push(tmp_row);
        }


        let n_r = data.len();
        let n_c = data[0].len();
        return Matrix { rows: n_r, cols: n_c, data };
    }

    pub fn copy(&self) -> Matrix {
        let mut n_data: Vec<Vec<f64>> = Vec::new();
        
        self.data.iter().for_each(|row| n_data.push(row.to_vec()));

        return Matrix { rows : self.rows, cols: self.cols, data: n_data };
    }

    pub fn print(&self) {
        self.data.iter().for_each(|v| println!("{:?}", v));
        println!();
    }

    pub fn identity(&mut self) {
        if self.rows != self.cols {
            panic!("Not a square matrix.");
        }
        for r in 0..self.rows {
            self.data[r][r] = 1.0;
        }
    }

    pub fn apply(&mut self, f: impl Fn(f64) -> f64) {
        self.data = self.data.iter()
                            .map(|v| {
                                v.iter()
                                .map(|x| f(*x))
                                .collect()
                            })
                            .collect();
    }

    pub fn combine(&self, b: Matrix, f: impl Fn(f64, f64) -> f64) -> Matrix {
        if self.rows != b.rows || self.cols != b.cols {
            panic!("Matrices must be of the same size");
        }
        let mut new_matrix = Matrix::new(self.rows, self.cols);
        for r in 0..self.rows {
            new_matrix.data[r] =
                self.data[r].iter().zip(b.data[r].iter()).map(|(a,b)| f(*a,*b)).collect();
        }
        return new_matrix;
    }


    pub fn dot(&self, b: Matrix) -> Matrix {
        if self.rows != b.cols || self.cols != b.rows {
            panic!("Dimensions not matched. M1 is {} by {}, M2 is {} by {}.", self.rows, self.cols, b.rows, b.cols);
        }
        let mut dp = Matrix::new(self.rows, b.cols);
        for i in 0..self.rows {
            for j in 0..b.cols {
                let mut sum = 0.0;
                for k in 0..b.rows {
                    sum += self.data[i][k] * b.data[k][j];
                }
                dp.data[i][j] = sum;
            }
        }
        return dp;
    }

    pub fn rref(&mut self) {
        if self.data[0][0] == 0.0 {
            swap_rows(self, 0);
        }
        let mut lead: usize = 0;
        let rows = self.rows;
        while lead < rows {
            for r in 0..rows {
                let div = self.data[lead][lead];
                let mult = self.data[r][lead] / div;

                if r == lead {
                    self.data[lead] = self.data[lead]
                                        .iter()
                                        .map(|entry| entry / div)
                                        .collect();
                }
                else {
                    for c in 0..self.cols {
                        self.data[r][c] -= self.data[lead][c] * mult;
                    }
                }
            }
            lead += 1;
        }
        correct(self);
    }

    pub fn cofactor(&self, expanded_row: usize, j: usize) -> f64 {
        let mut cut: Vec<Vec<f64>> = Vec::new();
        for r in 0..self.rows {
            if r == expanded_row {
                continue;
            }
            let mut v: Vec<f64> = Vec::new();
            for c in 0..self.cols {
                if c == j {
                    continue;
                }
                v.push(self.data[r][c]);
            }
            cut.push(v);
        }
        let n_r = cut.len();
        let n_c = cut[0].len();
        let minor = Matrix { rows : n_r, cols: n_c, data: cut }.det();
        let base: i32 = -1;
        return minor * f64::from(base.pow((expanded_row + j) as u32));
    }

    pub fn det(&self) -> f64 {
        if self.rows != self.cols {
            panic!("Determinant requires matrix to be a square. Input matrix was {:?}.", self);
        }
        if self.rows == 2 && self.cols == 2 {
            return self.data[0][0]*self.data[1][1] - self.data[0][1]*self.data[1][0];
        }
        else {
            let row: usize = 1;
            let mut det = 0.0;

            for j in 0..self.data[row].len() {
                det += self.cofactor(row, j) * self.data[row][j];
            }
            return det;
        }
    }

    pub fn transpose(&self) -> Matrix {
        let mut t = Matrix::new(self.cols, self.rows);
        for i in 0..self.rows {
            for j in 0..self.cols {
                t.data[j][i] = self.data[i][j];
            }
        }
        return t;
    }

    pub fn inverse(&self) -> Matrix {
        let d = self.det();
        if d == 0.0 {
            panic!("Determinant is 0. No inverse.");
        }

        let mut inv = Matrix::new(self.rows, self.cols);

        for row in 0..self.rows {
            for col in 0..self.cols {
                inv.data[row][col] = self.cofactor(row, col);
            }
        }

        correct(&mut inv);

        inv = inv.transpose();
        inv.apply(|x| x / d);
        return inv;
    }

}

fn swap_rows(m: &mut Matrix, row: usize){
    let mut n_r = 0;
    for r in 0..m.rows {
        if m.data[r][0] > 0.0 {
            n_r = r;
            break;
        }
    }
    let temp: Vec<f64> = m.data[row].clone();
    m.data[row] = m.data[n_r].clone();
    m.data[n_r] = temp;
}

fn correct(m: &mut Matrix) {
    for row in 0..m.rows {
        for col in 0..m.cols {
            let elem = m.data[row][col];
            if elem == -0.0 {
                m.data[row][col] = 0.0;
            }
            let floored = elem.floor();
            if elem - floored > 0.9999999 {
                m.data[row][col] = elem.round();
            }
            if elem > 0.0 && elem < 0.000001 {
                m.data[row][col] = 0.0;
            }
            if elem < 0.0 && elem > -0.00001 {
                m.data[row][col] = 0.0;
            }
        }
    }
}