use std::{fmt::Display, fs};

#[derive(Debug, PartialEq, Clone)]
pub struct Matrix {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<Vec<f64>>,
}

impl Matrix {
    pub fn new(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            data: vec![vec![0.0; cols]; rows],
        }
    }

    pub fn from_file(path: &str) -> Self {
        let content = fs::read_to_string(path).unwrap_or_else(|e| panic!("{e}"));
        let mut matrix: Vec<Vec<f64>> = Vec::new();
        for rows in content.lines() {
            let mut row: Vec<f64> = Vec::new();
            let entries: Vec<&str> = rows.split_whitespace().collect();

            entries
                .iter()
                .for_each(|ent| row.push(ent.parse::<f64>().unwrap()));

            matrix.push(row);
        }

        Self {
            rows: matrix.len(),
            cols: matrix[0].len(),
            data: matrix,
        }
    }

    pub fn from_string(input: &str) -> Self {
        let mut data: Vec<Vec<f64>> = Vec::new();
        let rows: Vec<&str> = input.split(';').collect();
        for row in rows {
            let entries: Vec<&str> = row.split_whitespace().collect();
            let mut tmp_row: Vec<f64> = Vec::new();

            entries
                .iter()
                .for_each(|ent| tmp_row.push(ent.parse::<f64>().unwrap()));

            data.push(tmp_row);
        }

        let n_r = data.len();
        let n_c = data[0].len();
        Self {
            rows: n_r,
            cols: n_c,
            data,
        }
    }

    pub fn copy(&self) -> Self {
        let mut n_data: Vec<Vec<f64>> = Vec::new();

        self.data.iter().for_each(|row| n_data.push(row.to_vec()));

        Self {
            rows: self.rows,
            cols: self.cols,
            data: n_data,
        }
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
        self.data = self
            .data
            .iter()
            .map(|v| v.iter().map(|x| f(*x)).collect())
            .collect();
    }

    pub fn combine(&self, b: Self, f: impl Fn(f64, f64) -> f64) -> Self {
        if self.rows != b.rows || self.cols != b.cols {
            panic!("Matrices must be of the same size");
        }
        let mut new_matrix = Self::new(self.rows, self.cols);
        for r in 0..self.rows {
            new_matrix.data[r] = self.data[r]
                .iter()
                .zip(b.data[r].iter())
                .map(|(a, b)| f(*a, *b))
                .collect();
        }
        new_matrix
    }

    pub fn dot(&self, b: Self) -> Self {
        if self.rows != b.cols || self.cols != b.rows {
            panic!(
                "Dimensions not matched. M1 is {} by {}, M2 is {} by {}.",
                self.rows, self.cols, b.rows, b.cols
            );
        }
        let mut dp = Self::new(self.rows, b.cols);
        for i in 0..self.rows {
            for j in 0..b.cols {
                let mut sum = 0.0;
                for k in 0..b.rows {
                    sum += self.data[i][k] * b.data[k][j];
                }
                dp.data[i][j] = sum;
            }
        }
        dp
    }

    pub fn rref(&mut self) {
        if self.data[0][0] == 0.0 {
            self.swap_rows(0);
        }
        let mut lead: usize = 0;
        let rows = self.rows;
        while lead < rows {
            for r in 0..rows {
                let div = self.data[lead][lead];
                let mult = self.data[r][lead] / div;

                if r == lead {
                    self.data[lead] = self.data[lead].iter().map(|entry| entry / div).collect();
                } else {
                    for c in 0..self.cols {
                        self.data[r][c] -= self.data[lead][c] * mult;
                    }
                }
            }
            lead += 1;
        }
        self.correct();
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
        let minor = Self {
            rows: n_r,
            cols: n_c,
            data: cut,
        }
        .det();
        let base: i32 = -1;
        minor * f64::from(base.pow((expanded_row + j) as u32))
    }

    pub fn det(&self) -> f64 {
        if self.rows != self.cols {
            panic!(
                "Determinant requires matrix to be a square. Input matrix was {:?}.",
                self
            );
        }
        if self.rows == 2 && self.cols == 2 {
            self.data[0][0] * self.data[1][1] - self.data[0][1] * self.data[1][0]
        } else {
            let row: usize = 1;
            let mut det = 0.0;

            for j in 0..self.data[row].len() {
                det += self.cofactor(row, j) * self.data[row][j];
            }
            det
        }
    }

    pub fn transpose(&self) -> Self {
        let mut t = Self::new(self.cols, self.rows);
        for i in 0..self.rows {
            for j in 0..self.cols {
                t.data[j][i] = self.data[i][j];
            }
        }
        t
    }

    pub fn inverse(&self) -> Self {
        let d = self.det();
        if d == 0.0 {
            panic!("Determinant is 0. No inverse.");
        }

        let mut inv = Self::new(self.rows, self.cols);

        for row in 0..self.rows {
            for col in 0..self.cols {
                inv.data[row][col] = self.cofactor(row, col);
            }
        }

        inv.correct();
        inv = inv.transpose();
        inv.apply(|x| x / d);
        inv
    }

    fn swap_rows(&mut self, row: usize) {
        let mut n_r = 0;
        for r in 0..self.rows {
            if self.data[r][0] > 0.0 {
                n_r = r;
                break;
            }
        }
        let temp: Vec<f64> = self.data[row].clone();
        self.data[row] = self.data[n_r].clone();
        self.data[n_r] = temp;
    }

    fn correct(&mut self) {
        for row in 0..self.rows {
            for col in 0..self.cols {
                let elem = self.data[row][col];
                if elem == -0.0 {
                    self.data[row][col] = 0.0;
                }
                let floored = elem.floor();
                if elem - floored > 0.9999999 {
                    self.data[row][col] = elem.round();
                }
                if elem > 0.0 && elem < 0.000001 {
                    self.data[row][col] = 0.0;
                }
                if elem < 0.0 && elem > -0.00001 {
                    self.data[row][col] = 0.0;
                }
            }
        }
    }
}

impl Display for Matrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for v in self.data.iter() {
            writeln!(f, "{:?}", v)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_string() {
        let m = Matrix::from_string("1 2 3 ; 4 5 6");
        let expected = Matrix {
            rows: 2,
            cols: 3,
            data: vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]],
        };

        assert!(m == expected);
    }

    #[test]
    fn test_display() {
        let m = Matrix::from_string("1 2 3 ; 4 5 6");
        
        assert_eq!("[1.0, 2.0, 3.0]\n[4.0, 5.0, 6.0]\n", m.to_string())
    }
}
