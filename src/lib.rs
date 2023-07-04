use std::ops::{Index, IndexMut};
use std::{fmt::Display, fs};

#[derive(Debug, PartialEq, Clone)]
pub struct Matrix {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<f64>,
}

impl Matrix {
    pub fn new(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            data: vec![0.0; rows * cols],
        }
    }

    pub fn from_file(path: &str) -> Self {
        let content = fs::read_to_string(path).unwrap_or_else(|e| panic!("{e}"));
        let mut data: Vec<f64> = Vec::new();
        let mut cols: usize = 0;
        let mut count: usize = 0;

        for r in content.lines() {
            let entries: Vec<&str> = r.split_whitespace().collect();
            let c = entries.len();
            if count > 0 && cols != c {
                panic!("Columns don't match");
            }
            cols = c;
            count += 1;

            let temp: Vec<f64> = entries
                .iter()
                .map(|ent| ent.parse::<f64>().unwrap())
                .collect();
            
            for item in temp {
                data.push(item);
            }
        }

        Self {
            rows: content.lines().collect::<Vec<_>>().len(),
            cols,
            data,
        }
    }

    pub fn from_string(input: &str) -> Self {
        let mut data: Vec<f64> = Vec::new();
        let rows: Vec<&str> = input.split(';').collect();
        let mut cols: usize = 0;
        let mut count: usize = 0;

        for r in &rows {
            let entries: Vec<&str> = r.split_whitespace().collect();
            let c = entries.len();
            if count > 0 && cols != c {
                panic!("Columns don't match");
            }
            cols = c;
            count += 1;

            let temp: Vec<f64> = entries
                .iter()
                .map(|ent| ent.parse::<f64>().unwrap())
                .collect();
            
            for item in temp {
                data.push(item);
            }
        }

        Self {
            rows: rows.len(),
            cols,
            data,
        }
    }

    pub fn copy(&self) -> Self {
        let mut n_data: Vec<f64> = Vec::new();

        self.data.iter().for_each(|elem| n_data.push(*elem));

        Self {
            rows: self.rows,
            cols: self.cols,
            data: n_data,
        }
    }

    pub fn print(&self) {
        for r in 0..self.rows {
            print!("[");
            for c in 0..self.cols {
                if c == self.cols - 1 { print!("{:.3}", self[r][c]); } else { print!("{:.3} ", self[r][c]); }
            }
            println!("]");
        }
    }

    pub fn identity(&mut self) {
        if self.rows != self.cols {
            panic!("Not a square matrix.");
        }
        for r in 0..self.rows {
            self[r][r] = 1.0;
        }
    }

    pub fn apply(&mut self, f: impl Fn(f64) -> f64) {
        self.data = self.data.iter().map(|elem| f(*elem)).collect()
    }

    pub fn combine(&self, b: Self, f: impl Fn(f64, f64) -> f64) -> Self {
        if self.rows != b.rows || self.cols != b.cols {
            panic!("Matrices must be of the same size.");
        }
        let mut new_matrix = Self::new(self.rows, self.cols);
        new_matrix.data = self.data
            .iter()
            .zip(b.data.iter())
            .map(|(a, b)| f(*a, *b))
            .collect();
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
                    sum += self[i][k] * b[k][j];
                }
                dp[i][j] = sum;
            }
        }
        dp
    }

    pub fn rref(&mut self) {
        if self[0][0] == 0.0 {
            self.swap_rows(0);
        }
        let mut lead: usize = 0;
        let rows = self.rows;
        while lead < rows {
            for r in 0..rows {
                let div = self[lead][lead];
                let mult = self[r][lead] / div;

                if r == lead {
                    // self[lead] = self[lead].iter().map(|entry| entry / div).collect::<Vec<_>>();
                    self[lead]
                        .iter_mut()
                        .for_each(|elem| *elem = (*elem) / div);
                } else {
                    for c in 0..self.cols {
                        self[r][c] -= self[lead][c] * mult;
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
                v.push(self[r][c]);
            }
            cut.push(v);
        }
        let flattened = cut.clone().into_iter().flatten().collect();
        let n_r = cut.len();
        let n_c = cut[0].len();
        let minor = Self {
            rows: n_r,
            cols: n_c,
            data: flattened,
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
            self[0][0] * self[1][1] - self[0][1] * self[1][0]
        } else {
            let row: usize = 1;
            let mut det = 0.0;

            for j in 0..self[row].len() {
                det += self.cofactor(row, j) * self[row][j];
            }
            det
        }
    }

    pub fn transpose(&self) -> Self {
        let mut t = Self::new(self.cols, self.rows);
        for i in 0..self.rows {
            for j in 0..self.cols {
                t[j][i] = self[i][j];
            }
        }
        t
    }

    pub fn trace(&self) -> f64 {
        if self.rows != self.cols {
            panic!(
                "Trace requires matrix to be square. Input matrix was {:?}.",
                self
            );
        }
        let mut t: f64 = 0.0;
        for i in 0..self.rows {
            t += self[i][i];
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
                inv[row][col] = self.cofactor(row, col);
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
            if self[r][0] > 0.0 {
                n_r = r;
                break;
            }
        }
        let temp: Vec<f64> = self[row].to_vec();
        for c in 0..self.cols {
            self[row][c] = self[n_r][c];
            self[n_r][c] = temp[n_r * self.cols + c];
        }
    }

    fn correct(&mut self) {
        for row in 0..self.rows {
            for col in 0..self.cols {
                let elem = self[row][col];
                if elem == -0.0 {
                    self[row][col] = 0.0;
                }
                let floored = elem.floor();
                if elem - floored > 0.9999999 {
                    self[row][col] = elem.round();
                }
                if elem > 0.0 && elem < 0.000001 {
                    self[row][col] = 0.0;
                }
                if elem < 0.0 && elem > -0.00001 {
                    self[row][col] = 0.0;
                }
            }
        }
    }
}

impl Display for Matrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for r in 0..self.rows {
            write!(f, "[")?;
            for c in 0..self.cols {
                if c == self.cols - 1 { write!(f, "{:.3}", self[r][c])?; } else { write!(f, "{:.3} ", self[r][c])?; }
                
            }
            writeln!(f, "]")?;
        }

        Ok(())
    }
}

impl Index<usize> for Matrix {
    type Output = [f64];

    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index * self.cols..(index + 1) * self.cols]
    }
}

impl IndexMut<usize> for Matrix {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index * self.cols..(index + 1) * self.cols]
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
            data: vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        };

        assert!(m == expected);
    }

    #[test]
    fn test_display() {
        let m = Matrix::from_string("1 2 3 ; 4 5 6");

        assert_eq!("[1.0, 2.0, 3.0]\n[4.0, 5.0, 6.0]\n", m.to_string())
    }
}
