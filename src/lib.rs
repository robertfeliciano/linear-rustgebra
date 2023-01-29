use std::{fmt::Display, fs};

#[derive(Debug, PartialEq, Clone)]
pub struct Matrix {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<f64>,
}
impl std::ops::Index<[usize; 2]> for Matrix {
    type Output = f64;
    fn index(&self, idx: [usize; 2]) -> &f64 {
        &self.data[idx[0] * self.cols + idx[1]]
    }
}
impl std::ops::IndexMut<[usize; 2]> for Matrix {
    fn index_mut(&mut self, idx: [usize; 2]) -> &mut f64 {
        &mut self.data[idx[0] * self.cols + idx[1]]
    }
}

impl Matrix {
    pub fn new(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            data: vec![0.0; rows * cols],
        }
    }
    fn parse_row(row: &str) -> Vec<f64> {
        row.split_whitespace()
            .map(|x| x.parse::<f64>().unwrap())
            .collect()
    }   
    
    pub fn from_string_by_sep(input: &str, f: impl Fn(&str) -> Vec<&str>) -> Self {
        let rows: Vec<&str> = f(input);
        let row_count = rows.len();
        let matrix: Vec<f64> = rows.iter().flat_map(|v| Self::parse_row(v)).collect();
        Self {
            rows: row_count,
            cols: matrix.len() / row_count,
            data: matrix,
        }
    }
    pub fn from_file(path: &str) -> Self {
        let content = fs::read_to_string(path).unwrap_or_else(|e| panic!("{e}"));
        Self::from_string_by_sep(&content, |x| x.lines().collect())        
    }
    pub fn from_string(input: &str) -> Self {
        Self::from_string_by_sep(input, |x| x.split(';').collect())        
    }

    pub fn copy(&self) -> Self {
        Self {
            rows: self.rows,
            cols: self.cols,
            data: self.data.clone(),
        }
    }
    pub fn print(&self) {
        self.data.chunks(self.cols).for_each(|x| println!("{x:?}"));
    }

    pub fn identity(&mut self) {
        if self.rows != self.cols {
            panic!("Not a square matrix.");
        }
        for r in 0..self.rows {
            self[[r, r]] = 1.0;
        }
    }

    pub fn apply(&mut self, f: impl Fn(f64) -> f64) {
        self.data = self.data.iter().map(|v| f(*v)).collect()
    }

    pub fn combine(&self, b: Self, f: impl Fn(f64, f64) -> f64) -> Self {
        if self.rows != b.rows || self.cols != b.cols {
            panic!("Matrices must be of the same size");
        }
        let data = self
            .data
            .iter()
            .zip(b.data.iter())
            .map(|(a, b)| f(*a, *b))
            .collect();
        Matrix {
            rows: self.rows,
            cols: self.cols,
            data,
        }
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
                    sum += self[[i, k]] * b[[i, k]];
                }
                dp[[i, j]] = sum;
            }
        }
        dp
    }
    pub fn rref(&mut self) {
        if self[[0, 0]] == 0.0 {
            self.swap_rows(0);
        }
        let mut lead: usize = 0;
        while lead < self.rows {
            for current_row in 0..self.rows {
                let div = self[[lead, lead]];
                let mult = self[[current_row, lead]] / div;
                for c in 0..self.cols {
                    if current_row == lead {
                        self[[lead, c]] /= div;
                    } else {
                        self[[current_row, c]] -= self[[lead, c]] * mult;
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
                v.push(self[[r, c]]);
            }
            cut.push(v);
        }
        let n_r = cut.len();
        let n_c = cut[0].len();

        let minor = Self {
            rows: n_r,
            cols: n_c,
            data: cut.into_iter().flatten().collect(),
        }
        .det();
        let base: i32 = -1;
        minor * f64::from(base.pow((expanded_row + j) as u32))
    }

    pub fn det(&self) -> f64 {
        if self.rows != self.cols {
            panic!("Determinant requires matrix to be a square. Input matrix was {self:?}.");
        }
        if self.rows == 2 && self.cols == 2 {
            self[[0, 0]] * self[[1, 1]] - self[[0, 1]] * self[[1, 0]]
        } else {
            let row: usize = 1;
            let mut det = 0.0;

            for j in 0..self.cols {
                det += self.cofactor(row, j) * self[[row, j]];
            }
            det
        }
    }

    pub fn transpose(&self) -> Self {
        let mut t = Self::new(self.cols, self.rows);
        for i in 0..self.rows {
            for j in 0..self.cols {
                t[[j, i]] = self[[i, j]];
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
                inv[[row, col]] = self.cofactor(row, col);
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
            if self[[r, 0]] > 0.0 {
                n_r = r;
                break;
            }
        }
        self.print();
        let (slice, slice2) = self.data.split_at_mut(row * self.rows + self.cols);
        let slice = &mut slice[row * self.rows..];
        let slice2 = &mut slice2[n_r * (self.rows - row)..n_r * (self.rows - row) + self.cols];
        slice.swap_with_slice(slice2);
        self.print();
    }

    fn correct(&mut self) {
        for row in 0..self.rows {
            for col in 0..self.cols {
                let elem = self[[row, col]];
                if elem == -0.0 {
                    self[[row, col]] = 0.0;
                }
                let floored = elem.floor();
                if elem - floored > 0.9999999 {
                    self[[row, col]] = elem.round();
                }
                if elem > 0.0 && elem < 0.000001 {
                    self[[row, col]] = 0.0;
                }
                if elem < 0.0 && elem > -0.00001 {
                    self[[row, col]] = 0.0;
                }
            }
        }
    }
}

impl Display for Matrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for v in self.data.iter() {
            writeln!(f, "{v:?}")?;
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
