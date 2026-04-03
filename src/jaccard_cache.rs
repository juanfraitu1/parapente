//! Jaccard computation cache for avoiding redundant calculations
//!
//! Uses an LRU cache to store pairwise Jaccard results.
//! Key insight: Jaccard(A,B) = Jaccard(B,A), so we canonicalize keys.

use hashbrown::HashMap;
use roaring::RoaringBitmap;
use std::sync::{Arc, Mutex};

/// Cache key for Jaccard computations (canonical ordering)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct JaccardKey(u64, u64);

impl JaccardKey {
    /// Create a canonical key from two bitmap IDs (order doesn't matter)
    fn new(a: u64, b: u64) -> Self {
        if a <= b {
            Self(a, b)
        } else {
            Self(b, a)
        }
    }
}

/// Thread-safe Jaccard computation cache
pub struct JaccardCache {
    cache: Mutex<HashMap<JaccardKey, f64>>,
    hits: Mutex<u64>,
    misses: Mutex<u64>,
}

impl JaccardCache {
    /// Create a new empty cache
    pub fn new() -> Self {
        Self {
            cache: Mutex::new(HashMap::with_capacity(10000)),
            hits: Mutex::new(0),
            misses: Mutex::new(0),
        }
    }
    
    /// Get or compute Jaccard similarity between two bitmaps
    pub fn get_or_compute(&self, id1: u64, bm1: &RoaringBitmap, id2: u64, bm2: &RoaringBitmap) -> f64 {
        let key = JaccardKey::new(id1, id2);
        
        // Fast path: check cache
        {
            let cache = self.cache.lock().unwrap();
            if let Some(&value) = cache.get(&key) {
                *self.hits.lock().unwrap() += 1;
                return value;
            }
        }
        
        // Compute Jaccard
        *self.misses.lock().unwrap() += 1;
        let jaccard = Self::compute_jaccard(bm1, bm2);
        
        // Store in cache
        {
            let mut cache = self.cache.lock().unwrap();
            cache.insert(key, jaccard);
        }
        
        jaccard
    }
    
    /// Compute Jaccard directly without caching
    #[inline]
    pub fn compute_jaccard(bm1: &RoaringBitmap, bm2: &RoaringBitmap) -> f64 {
        if bm1.is_empty() && bm2.is_empty() {
            return 0.0;
        }
        
        let intersection = bm1.intersection_len(bm2);
        if intersection == 0 {
            return 0.0;
        }
        
        let union = bm1.union_len(bm2);
        if union == 0 {
            0.0
        } else {
            intersection as f64 / union as f64
        }
    }
    
    /// Get cache statistics
    pub fn stats(&self) -> (u64, u64, f64) {
        let hits = *self.hits.lock().unwrap();
        let misses = *self.misses.lock().unwrap();
        let total = hits + misses;
        let hit_rate = if total > 0 {
            hits as f64 / total as f64
        } else {
            0.0
        };
        (hits, misses, hit_rate)
    }
    
    /// Clear the cache
    pub fn clear(&self) {
        self.cache.lock().unwrap().clear();
        *self.hits.lock().unwrap() = 0;
        *self.misses.lock().unwrap() = 0;
    }
    
    /// Get cache size
    pub fn len(&self) -> usize {
        self.cache.lock().unwrap().len()
    }
}

impl Default for JaccardCache {
    fn default() -> Self {
        Self::new()
    }
}

/// Global Jaccard cache instance
lazy_static::lazy_static! {
    pub static ref GLOBAL_JACCARD_CACHE: Arc<JaccardCache> = Arc::new(JaccardCache::new());
}

/// Cached Jaccard computation using global cache
#[inline]
pub fn cached_jaccard(id1: u64, bm1: &RoaringBitmap, id2: u64, bm2: &RoaringBitmap) -> f64 {
    GLOBAL_JACCARD_CACHE.get_or_compute(id1, bm1, id2, bm2)
}

/// Report global cache statistics
pub fn report_cache_stats() {
    let (hits, misses, hit_rate) = GLOBAL_JACCARD_CACHE.stats();
    let total = hits + misses;
    if total > 0 {
        println!("  Jaccard cache: {} hits, {} misses ({:.1}% hit rate)",
                 hits, misses, hit_rate * 100.0);
    }
}
