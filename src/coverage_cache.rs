//! LRU cache for coverage profiles to avoid recomputation
//!
//! Many loci overlap or are queried multiple times - cache the coverage profiles.

use lru::LruCache;
use std::num::NonZeroUsize;
use std::sync::Mutex;

/// Cache key: (chrom, start, end, bin_size)
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct CoverageKey(String, i64, i64, i64);

/// Cached coverage profile
#[derive(Debug, Clone)]
pub struct CachedCoverage {
    pub bins: Vec<f64>,
    pub max: f64,
    pub mean: f64,
    pub std: f64,
}

/// Thread-safe LRU cache for coverage profiles
pub struct CoverageCache {
    cache: Mutex<LruCache<CoverageKey, CachedCoverage>>,
    hits: Mutex<u64>,
    misses: Mutex<u64>,
}

impl CoverageCache {
    /// Create a new cache with the given capacity
    pub fn new(capacity: usize) -> Self {
        let capacity = NonZeroUsize::new(capacity).unwrap_or(NonZeroUsize::new(1000).unwrap());
        Self {
            cache: Mutex::new(LruCache::new(capacity)),
            hits: Mutex::new(0),
            misses: Mutex::new(0),
        }
    }
    
    /// Get a cached coverage profile if it exists
    pub fn get(&self, chrom: &str, start: i64, end: i64, bin_size: i64) -> Option<CachedCoverage> {
        let key = CoverageKey(chrom.to_string(), start, end, bin_size);
        let result = self.cache.lock().unwrap().get(&key).cloned();
        
        if result.is_some() {
            *self.hits.lock().unwrap() += 1;
        }
        
        result
    }
    
    /// Insert a coverage profile into the cache
    pub fn insert(&self, chrom: &str, start: i64, end: i64, bin_size: i64, coverage: CachedCoverage) {
        let key = CoverageKey(chrom.to_string(), start, end, bin_size);
        *self.misses.lock().unwrap() += 1;
        self.cache.lock().unwrap().put(key, coverage);
    }
    
    /// Get or compute coverage profile
    pub fn get_or_compute<F>(&self, chrom: &str, start: i64, end: i64, bin_size: i64, compute: F) -> CachedCoverage
    where
        F: FnOnce() -> anyhow::Result<(Vec<f64>, f64, f64, f64)>,
    {
        // Try cache first
        if let Some(cached) = self.get(chrom, start, end, bin_size) {
            return cached;
        }
        
        // Compute
        match compute() {
            Ok((bins, max, mean, std)) => {
                let coverage = CachedCoverage { bins, max, mean, std };
                self.insert(chrom, start, end, bin_size, coverage.clone());
                coverage
            }
            Err(_) => {
                // Return empty on error
                CachedCoverage {
                    bins: vec![],
                    max: 0.0,
                    mean: 0.0,
                    std: 0.0,
                }
            }
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
    
    /// Get current cache size
    pub fn len(&self) -> usize {
        self.cache.lock().unwrap().len()
    }
}

impl Default for CoverageCache {
    fn default() -> Self {
        Self::new(2000)
    }
}

use std::sync::Arc;

/// Global coverage cache instance
lazy_static::lazy_static! {
    pub static ref GLOBAL_COVERAGE_CACHE: Arc<CoverageCache> = Arc::new(CoverageCache::default());
}

/// Report global cache statistics
pub fn report_coverage_cache_stats() {
    let (hits, misses, hit_rate) = GLOBAL_COVERAGE_CACHE.stats();
    let total = hits + misses;
    if total > 0 {
        println!("  Coverage cache: {} hits, {} misses ({:.1}% hit rate)",
                 hits, misses, hit_rate * 100.0);
    }
}
