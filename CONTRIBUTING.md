# Contributing to reciprocal_spouge

Thank you for your interest in this project. Contributions are most welcome! Whether it is a bug fix, new feature, or improvements to documentation.

## Getting Started
1. **Fork the repository and then clone it repository to your local machine** `git clone https://github.com`
2. **Install the Rust toolchain** using [rustup](http://rustup.rs). We use the stable toolchain: Run `rust update stable`
3. **Create a new branch for your work** `git checkout -b features/my-features`.

## Development Workflow
Before making a submission, please ensure your code aligns with our projects standards:
* **Formatting:** Run `cargo fmt` to auto-format your code.
* **Linting:** Run `cargo clippy --all-targets --all-features` and resolve all Clippy errors (and often warnings) that trigger your Continuous Integration (CI) pipeline to fail.
* **Testing:** Run `cargo test` to verify all tests pass.

## How to Suggest New Features
1. **Search Existing Issues** Check if someone else has already suggested the feature.
2. **Open a Feature Request** Use our "Feature Request" issue template to describe:
     * The problem the proposed feature solves.
     * Your proposed solution or API design.
     * Alternative workarounds you have considered.
3. **Wait for Feedback** Let the maintainers discuss the design with you. This prevents time being wasted on code that might not necessarily align with the project's goals.

## Documentation
1. Update your branch to the latest `main` with `git pull origin main`.
2. Push your changes to your fork and [create a Pull Request on Github]
3. Ensure your PR description clearly states the problem solved and links relevant issues.

## Code of Conduct
We abide by the Contributor Convenant Code of Conduct [Link to CODE_OF_CONDUCT.md]. By participating in this project you agree to do the same. Please be respectful and collaborative.

## Getting Help
TODO
