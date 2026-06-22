//! Build libhex9's static library (CMake `hex9Static` target) and generate
//! Rust FFI bindings from the C ABI header `hex9_c.h`.
//!
//! The library is linked statically and carries no external runtime dependency
//! (the static target deliberately omits OpenMP; the embedded authalic-warp
//! blob is `.incbin`'d at compile time). The addressing layout follows libhex9's
//! CMake default — L30 (the reclaimed layout) unless built with HEX9_USE_L29=ON.

use std::env;
use std::path::PathBuf;

fn main() {
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    // The libhex9 repository root is this crate's parent directory.
    let libhex9_root = manifest
        .parent()
        .expect("hex9-sys must live inside the libhex9 repo (libhex9/geoplegma/)")
        .to_path_buf();
    let header = libhex9_root.join("hex9_c.h");

    // 1. Build the static library via CMake. HEX9_PYTHON=OFF skips the nanobind
    //    module; build_target builds only the archive (and its deps), not tests.
    let dst = cmake::Config::new(&libhex9_root)
        .define("HEX9_PYTHON", "OFF")
        .build_target("hex9Static")
        .build();
    let build_dir = dst.join("build");

    println!("cargo:rustc-link-search=native={}", build_dir.display());
    // Some generators nest per-config; add the common Release subdir too.
    println!(
        "cargo:rustc-link-search=native={}",
        build_dir.join("Release").display()
    );
    println!("cargo:rustc-link-lib=static=hex9");

    // libhex9's core is C++, so pull in the C++ runtime.
    let target_os = env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();
    if target_os == "macos" || target_os == "ios" {
        println!("cargo:rustc-link-lib=dylib=c++");
    } else {
        println!("cargo:rustc-link-lib=dylib=stdc++");
    }

    // 2. Generate bindings from the C ABI surface.
    let bindings = bindgen::Builder::default()
        .header(header.to_str().unwrap())
        .allowlist_function("hex9_.*")
        .allowlist_type("hex9_.*")
        .allowlist_var("HEX9_.*")
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .generate()
        .expect("failed to generate libhex9 FFI bindings");

    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_dir.join("bindings.rs"))
        .expect("failed to write libhex9 bindings");

    println!("cargo:rerun-if-changed={}", header.display());
}
