//! Intern read names to dense u32 ids for fast set operations and smaller indices.

use fxhash::FxHashSet as HashSet;
use hashbrown::HashMap;
use roaring::RoaringBitmap;

/// Maps read name strings to stable u32 ids and back for decoding at API boundaries.
#[derive(Debug, Default)]
pub struct ReadInterner {
    forward: HashMap<String, u32>,
    backward: Vec<String>,
}

impl ReadInterner {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn intern(&mut self, name: &str) -> u32 {
        if let Some(&id) = self.forward.get(name) {
            return id;
        }
        let id = self.backward.len() as u32;
        self.backward.push(name.to_string());
        self.forward.insert(name.to_string(), id);
        id
    }

    pub fn id_to_string(&self, id: u32) -> Option<&str> {
        self.backward.get(id as usize).map(|s| s.as_str())
    }

    /// Number of interned read names.
    pub fn len(&self) -> usize {
        self.backward.len()
    }

    pub fn bitmap_to_hashset(&self, bm: &RoaringBitmap) -> HashSet<String> {
        let mut out = HashSet::default();
        out.reserve(bm.len() as usize);
        for id in bm.iter() {
            if let Some(s) = self.backward.get(id as usize) {
                out.insert(s.clone());
            }
        }
        out
    }
}
