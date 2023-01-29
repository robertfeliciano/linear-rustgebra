use linalg::Matrix;

fn main() {
    let m1 = Matrix::from_string("1 2 3 ; 4 5 6");
    println!("{m1}");
    let m2 = Matrix::from_string("7 8 9; 10 11 12");
    println!("{m2}");
    let m3 = m1.combine(m2, |a, b| (a * b));
    println!("{m3}");

    println!("{}",m3[[5,5]]);

    // let mut m: Matrix = Matrix::new(3,3);
    // m.identity();

    // let mut mcpy = m.copy();
    // mcpy.apply(|x| x+3.0);
    // mcpy.print();

    // m.print();

    // let m1 = Matrix::from_file("src/2b2.txt");
    // m1.print();
    // println!("det(m1) = {}", m1.det());

    // let mut m2 = Matrix::from_file("src/m1.txt");
    // m2.apply(|i| i*2.0);
    // m2.print();
    // println!("det(m2) = {}", m2.det());

    // let m3 = Matrix::from_string("1 2 3 ; 4 5 6 ; 7 8 9");
    // m3.print();
    // println!("det(m3) = {}", m3.det());

    // let m4 = Matrix::from_string("9 8 4 ; 2 3 7; 4 1 1");
    // m4.print();
    // println!("det(m4) = {}", m4.det());

    // let m34 = m4.dot(m3);
    // println!("m4 dot m3:");
    // m34.print();

    // let mut m4t = m4.transpose();
    // println!("testing apply");
    // m4t.print();
    // m4t.apply(|x| x+99.0);
    // m4t.print();

    // let m5 = Matrix::from_string("3 0 2 ; 2 0 -2 ; 0 1 1 ");
    // m5.print();

    // let m5i = m5.inverse();
    // println!("inverse of m5:");
    // m5i.print();

    // // let mut m6 = Matrix::from_string("0 1 2 ; 0 3 1 ; 5 2 2");
    // let mut m6 = Matrix::from_string("5 -6 -7 7 ;
    //                                                 3 -2 5 -17 ;
    //                                                 2 4 -3 29");
    // m6.print();
    // println!("Row Reduce Echelon Form calculation:");
    // // println!("{:?}", m6);
    // m6.rref();
    // m6.print();

    // let mut m7 = Matrix::from_file("src/4b5.txt");
    // m7.print();
    // println!("{:?}", m7);
    // println!("Row Reduce Echelon Form calculation:");
    // m7.rref();
    // m7.print();
}
