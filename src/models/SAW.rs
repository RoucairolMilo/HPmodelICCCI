extern crate nalgebra;
use std::collections::HashMap;
use nalgebra::{ArrayStorage, Const, Matrix, Vector3};

#[derive(PartialEq, Eq, Hash, Clone, Copy)]
pub struct Move{
    pub dir : Vector3<i32>,
    pub dir2 : Vector3<i32>,
    pub dir3 : Vector3<i32>,
    //pub len : usize
}

const up: Matrix<i32, Const<3>, Const<1>, ArrayStorage<i32, 3, 1>> = Vector3::new(0, 1, 0);
const down: Matrix<i32, Const<3>, Const<1>, ArrayStorage<i32, 3, 1>> = Vector3::new(0, -1, 0);
const left: Matrix<i32, Const<3>, Const<1>, ArrayStorage<i32, 3, 1>> = Vector3::new(-1, 0, 0);
const right: Matrix<i32, Const<3>, Const<1>, ArrayStorage<i32, 3, 1>> = Vector3::new(1, 0, 0);
const fwd: Matrix<i32, Const<3>, Const<1>, ArrayStorage<i32, 3, 1>> = Vector3::new(0, 0, 1);
const bwd: Matrix<i32, Const<3>, Const<1>, ArrayStorage<i32, 3, 1>> = Vector3::new(0, 0, -1);
const sixDir: [Matrix<i32, Const<3>, Const<1>, ArrayStorage<i32, 3, 1>>; 6] = [up, down, left, right, fwd, bwd];

const sizeX : i32 = 10;
const sizeY : i32 = 10;
const sizeZ : i32 = 10;

#[derive(Clone, Copy)]
pub struct Cell{
    pub prev : Vector3<i32>, //début si c'est 0, 0, 0
    pub next : Vector3<i32>  //fin si c'est 1000, 1000, 1000
}

#[derive(Clone)]
pub struct State{
    pub lattice : HashMap<Vector3<i32>, Cell>,
    pub anchor : Vector3<i32>,
    pub seq : Vec<Move>,
}

impl State{
    pub const CONSIDER_NON_TERM: bool = false;

    pub fn new() -> Self {
        Self{
            lattice: HashMap::new(),
            anchor: Vector3::new(0, 0, 0),
            seq : Vec::new()
        }
    }

    pub fn play(&mut self, m : Move){
        let pt = self.anchor+m.dir;
        if !self.lattice.contains_key(&pt) {
            let c = Cell{
                prev: self.anchor.clone(),
                next: Vector3::new(1000, 1000, 1000)
            };

            let mut parent = self.getFromLattice(self.anchor).clone();
            parent.next = pt;
            self.lattice.insert(self.anchor, parent);
            self.lattice.insert(pt, c);
            self.anchor = pt;
            self.seq.push(m);
        }
    }

    pub fn legal_moves(& self) -> Vec<Move>{
        let mut vec :Vec<Move> = Vec::new();
        let mut dir2 = Vector3::new(0, 0, 0);
        let mut dir3 = Vector3::new(0, 0, 0);

        if self.seq.len() >= 1 {
            dir2 = self.seq[self.seq.len()-1].dir;
            dir3 = self.seq[self.seq.len()-1].dir2;
        }

        for dir in sixDir {
            let pt = self.anchor+dir;
            if !self.lattice.contains_key(&pt) && self.inLattice(pt){
                let mv = Move{dir : dir, dir2 : dir2, dir3 : dir3 };
                //let mv = Move{dir : dir, len : self.seq.len()};
                vec.push(mv);
            }
        }
        return vec;
    }

    pub fn inLattice(&self, pt : Vector3<i32>) -> bool{
        return  pt[0] >=0 && pt[0] < sizeX && pt[1] >=0 && pt[1] < sizeY && pt[2] >=0 && pt[2] < sizeZ;
    }

    pub fn terminal(& self) -> bool{
        return self.legal_moves().len() == 0;
    }

    pub fn score(& self) -> f64 {

        return self.seq.len() as f64;
    }

    pub fn heuristic(& self, m : Move) -> f64{
        return 0.0;
    }

    pub fn smoothedScore(&self) ->f64{
        return self.score();
    }

    pub fn getFromLattice(& self, key : Vector3<i32>) -> Cell {
        match self.lattice.get(&key) {
            Some(&cell) => return cell,
            None => return Cell{
                prev: Default::default(),
                next: Default::default()
            }.clone()
        }
    }

    pub fn printLattice(& self){
        for i in self.lattice.keys() {
            println!("{}", i);
        }
    }
}