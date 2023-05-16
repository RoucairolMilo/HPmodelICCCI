
extern crate nalgebra;
use std::collections::HashMap;
use nalgebra::{ArrayStorage, Const, Matrix, Vector3};

#[derive(PartialEq, Eq, Hash, Clone, Copy)]
pub struct Move{
    pub dir : Vector3<i32>,
    //pub dir2 : Vector3<i32>,
    //pub dir3 : Vector3<i32>,
    pub len : usize
}

const up: Matrix<i32, Const<3>, Const<1>, ArrayStorage<i32, 3, 1>> = Vector3::new(0, 1, 0);
const down: Matrix<i32, Const<3>, Const<1>, ArrayStorage<i32, 3, 1>> = Vector3::new(0, -1, 0);
const left: Matrix<i32, Const<3>, Const<1>, ArrayStorage<i32, 3, 1>> = Vector3::new(-1, 0, 0);
const right: Matrix<i32, Const<3>, Const<1>, ArrayStorage<i32, 3, 1>> = Vector3::new(1, 0, 0);
const fwd: Matrix<i32, Const<3>, Const<1>, ArrayStorage<i32, 3, 1>> = Vector3::new(0, 0, 1);
const bwd: Matrix<i32, Const<3>, Const<1>, ArrayStorage<i32, 3, 1>> = Vector3::new(0, 0, -1);
const sixDir: [Matrix<i32, Const<3>, Const<1>, ArrayStorage<i32, 3, 1>>; 6] = [up, down, left, right, fwd, bwd];

//1.0, -0.2, 0.0 de base
const HHrew : f64 = 1.0;
const HPrew : f64 = -0.2; //normalement c'est 0 mais ça aide pour le potentiel local de le mettre bas, pour éviter les liaisons HP car elles bloquent du potentiel
const PPrew : f64 = 0.0;

#[derive(Clone, Copy)]
pub struct Cell{
    pub monomere : char,
    pub prev : Vector3<i32>, //début si c'est 0, 0, 0
    pub next : Vector3<i32>  //fin si c'est 1000, 1000, 1000
}

#[derive(Clone)]
pub struct State{
    pub lattice : HashMap<Vector3<i32>, Cell>,
    pub anchor : Vector3<i32>,
    pub molecule : String,
    pub best_possible_score : f64,
    pub seq : Vec<Move>,
    pub reached_best_score :bool
}

pub struct mvl_st{ //pour iterative depth greedy moves
    pub mvl : Vec<Move>,
    pub st : State
}

/*

1 : HPHHPPHHHHPHHHPPHHPPHPHHHPHPHHPPHHPPPHPPPPPPPPHH - 32
2 : HHHHPHHPHHHHHPPHPPHHPPHPPPPPPHPPHPPPHPPHHPPHHHPH - 34
3 : PHPHHPHHHHHHPPHPHPPHPHHPHPHPPPHPPHHPPHHPPHPHPPHP - 34
4 : PHPHHPPHPHHHPPHHPHHPPPHHHHHPPHPHHPHPHPPPPHPPHPHP - 33
5 : PPHPPPHPHHHHPPHHHHPHHPHHHPPHPHPHPPHPPPPPPHHPHHPH - 32
6 : HHHPPPHHPHPHHPHHPHHPHPPPPPPPHPHPPHPPPHPPHHHHHHPH - 32
7 : PHPPPPHPHHHPHPHHHHPHHPHHPPPHPHPPPHHHPPHHPPHHPPPH - 32
8 : PHHPHHHPHHHHPPHHHPPPPPPHPHHPPHHPHPPPHHPHPHPHHPPP - 31
9 : PHPHPPPPHPHPHPPHPHHHHHHPPHHHPHPPHPHHPPHPHHHPPPPH - 34
10 : PHHPPPPPPHHPPPHHHPHPPHPHHPPHPPHPPHHPPHHHHHHHPPHH - 33


S1 : HPHPPHHPHPPHPHHPPHPH - 11
S2 : HHPPHPPHPPHPPHPPHPPHPPHH - 13
S3 : PPHPPHHPPPPHHPPPPHHPPPPHH - 9
S4 : PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP - 18
S5 : PPHPPHHPPHHPPPPPHHHHHHHHHHPPPPPPHHPPHHPPHPPHHHHH - 31
S6 : HHPHPHPHPHHHHPHPPPHPPPHPHHHHPHPHPHPHH - 31
S6 : HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH - 31
S7 : PPHHHPHHHHHHHHPPPHHHHHHHHHHPHPPPHHHHHHHHHHHHPPPPHHHHHHPHHPHP - 52
S8 : HHHHHHHHHHHHPHPHPPHHPPHHPPHPPHHPPHHPPHPPHHPPHHPPHPHPHHHHHHHHHHHH - 56

S8 :
best score yet : 51 after 1.7459986
best score yet : 54 after 1.7787958000000001
best score yet : 55 after 64.4129103
best score yet : 57 after 190.7638882
best score yet : 58 after 190.8340144


 */

impl State{
    pub const CONSIDER_NON_TERM: bool = false;

    /*
    pub fn new() -> Self {
        let mut mol = "HPHHPPHHHHPHHHPPHHPPHPHHHPHPHHPPHHPPPHPPPPPPPPHH".to_string();
        //"HHHPPPHHPHPHHPHHPHHPHPPPPPPPHPHPPHPPPHPPHHHHHHPH".chars().rev().collect::<String>();
        Self{
            lattice: HashMap::new(),
            anchor: Vector3::new(0, 0, 0),
            molecule : mol,
            best_possible_score : 32.0,
            seq : Vec::new(),
            reached_best_score : false
        }
    }

     */


    pub fn new() -> Self {
        let mut mol = "PPHPPHHPPHHPPPPPHHHHHHHHHHPPPPPPHHPPHHPPHPPHHHHH".to_string();
        //"HHHPPPHHPHPHHPHHPHHPHPPPPPPPHPHPPHPPPHPPHHHHHHPH".chars().rev().collect::<String>();
        let mut s = Self{
            lattice: HashMap::new(),
            anchor: Vector3::new(0, 0, 0),
            molecule : mol,
            best_possible_score : 140.0,
            seq : Vec::new(),
            reached_best_score : false
        };
        //s.playNoSeqPush(s.legal_moves()[0]);
        //s.playNoSeqPush(s.legal_moves()[0]);
        return s
    }

    pub fn playNoSeqPush(&mut self, m : Move){
        let pt = self.anchor+m.dir;
        if !self.lattice.contains_key(&pt) {
            let monomere = self.molecule.chars().nth(0).unwrap();
            let c = Cell{
                monomere: monomere,
                prev: self.anchor.clone(),
                next: Vector3::new(1000, 1000, 1000)
            };

            let mut parent = self.getFromLattice(self.anchor).clone();
            parent.next = pt;
            self.lattice.insert(self.anchor, parent);

            self.lattice.insert(pt, c);
            self.anchor = pt;
            self.molecule.remove(0);
            /*
            if monomere == 'H'{
                for d in sixDir{
                    if d != -m.dir && self.getFromLattice(pt + d).monomere  == 'H'{
                        self.score += HHrew;
                    }
                }
            }
            */

        }
    }


    pub fn play(&mut self, m : Move){
        let pt = self.anchor+m.dir;
        if !self.lattice.contains_key(&pt) {
            let monomere = self.molecule.chars().nth(0).unwrap();
            let c = Cell{
                monomere: monomere,
                prev: self.anchor.clone(),
                next: Vector3::new(1000, 1000, 1000)
            };

            let mut parent = self.getFromLattice(self.anchor).clone();
            parent.next = pt;
            self.lattice.insert(self.anchor, parent);

            self.lattice.insert(pt, c);
            self.anchor = pt;
            self.molecule.remove(0);
            self.seq.push(m);
            /*
            if monomere == 'H'{
                for d in sixDir{
                    if d != -m.dir && self.getFromLattice(pt + d).monomere  == 'H'{
                        self.score += HHrew;
                    }
                }
            }
            */

        }
    }

    pub fn legal_moves(& self) -> Vec<Move>{
        let mut vec :Vec<Move> = Vec::new();
        if self.molecule.len() == 0 {
            return vec;
        }

        let mut dir2 = Vector3::new(0, 0, 0);
        let mut dir3 = Vector3::new(0, 0, 0);

        if self.seq.len() >= 1 {
            //dir2 = self.seq[self.seq.len()-1].dir;
            //dir3 = self.seq[self.seq.len()-1].dir2;
        }

        for dir in sixDir {
            let pt = self.anchor+dir;
            if !self.lattice.contains_key(&pt) {
                //let mv = Move{dir : dir, dir2 : dir2, dir3 : dir3 };
                let mv = Move{dir : dir, len : self.molecule.len()};
                vec.push(mv);
            }
        }
        return vec;
    }


    pub fn iterative_depth_greedy_moves(&mut self) -> Vec<Move>{


        let mut open: Vec<mvl_st> = Vec::new();

        open.push(mvl_st{mvl : Vec::new(), st : self.clone()});
        //println!("hey {} {}", self.seq.len(), self.score());
        while open.len() != 0 {
            //println!("etape");
            let mut new_open: Vec<mvl_st>= Vec::new();
            for mut o in open {
                for m in o.st.legal_moves() {


                    if o.st.localPotential(o.st.anchor, m, o.st.molecule.chars().nth(0).unwrap()) == 1.0 {
                        o.mvl.push(m);
                        //println!("quoiiii");

                        return o.mvl
                    }

                    let mut mvl = o.mvl.clone();
                    mvl.push(m);
                    let mut st = o.st.clone();
                    st.play(m);
                    new_open.push(mvl_st{mvl : mvl, st : st})

                }
            }
            open = new_open;
        }

        //println!("hoho bug ?");
        let mut vec :Vec<Move> = Vec::new();
        vec.push(self.legal_moves()[0]);
        return vec;
    }

    pub fn terminal(& self) -> bool{
        return self.molecule == "";
    }

    pub fn score(&mut self) -> f64 {

        //return self.score;

        let mut score = 0.0;
        for i in self.lattice.keys() {
            for j in self.lattice.keys() {
                if manhattanDist(*i, *j) == 1{
                    let temp = self.getFromLattice(*i);
                    let temp2 = self.getFromLattice(*j);
                    if temp.prev != *j && temp.next != *j {
                        if temp.monomere == temp2.monomere {
                            if temp.monomere == 'H'{
                                score += HHrew;
                            }
                        }
                    }
                }
            }
        }
        if score/2.0 >= self.best_possible_score { self.reached_best_score = true }
        return score/2.0;
    }

    pub fn localPotential(& self, pos : Vector3<i32>, mov : Move, monomere : char) -> f64{

        let posmov = pos + mov.dir;
        if self.lattice.contains_key(&posmov){
            return -1.0;
        }

        let mut sum = 0.0;
        for dir in sixDir {
            let pt = posmov + dir;
            if self.lattice.contains_key(&pt) && pt != pos{
                if self.getFromLattice(pt).monomere != monomere{
                    sum += HPrew;
                }else{
                    if monomere == 'H' {
                        sum += HHrew;
                    }else{
                        sum += PPrew;
                    }
                }
            }
        }
        return sum;
    }

    pub fn heuristic(& self, m : Move) -> f64{
        if self.molecule.len() == 0 {
            //println!("ne devrait pas arriver");
            //return 0.0;
        }

        return self.localPotential(self.anchor, m, self.molecule.chars().nth(0).unwrap())
    }

    pub fn smoothedScore(&self) ->f64{
        let mut score = 0.0;
        for i in self.lattice.keys() {
            for j in self.lattice.keys() {
                if manhattanDist(*i, *j) == 1{
                    let temp = self.getFromLattice(*i);
                    let temp2 = self.getFromLattice(*j);
                    if temp.prev != *j && temp.next != *j {
                        if temp.monomere == temp2.monomere {
                            if temp.monomere == 'H'{
                                score += HHrew;
                            }else{
                                score += PPrew;
                            }
                        }else{
                            score += HPrew;
                        }
                    }
                }
            }
        }

        if score < 0.0 {return 0.0}
        return score/2.0;}

    pub fn getFromLattice(& self, key : Vector3<i32>) -> Cell {
        match self.lattice.get(&key) {
            Some(&cell) => return cell,
            None => return Cell{
                monomere: 'A',
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

pub fn manhattanDist(from : Vector3<i32>, to : Vector3<i32>) -> i32{
    let val = (from.x - to.x).abs() + (from.y - to.y).abs() + (from.z - to.z).abs();
    return val;
}