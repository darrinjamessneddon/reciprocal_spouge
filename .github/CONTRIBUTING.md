**CONTRIBUTING GUIDELINES**

Thank you for contributing to reciprocal_spouge! We welcome contributions from developers of all experience levels- whether it's a bug fix, testing, a new feature, or improvements to documentation.

**Code of Conduct**

This project abides by the Contributor Covenant Code of Conduct.
By participating in this project, you agree to abide by the terms of the `[CODE_OF_CONDUCT.md]`. Please be respectful and collaborative.

**How to Suggest New Features**

1. Search Existing Issues: Check if the feature has already been suggested by someone else.
2. Open a Feature Request: Use our "Feature Request" issue template to describe:
    * The problem this feature solves.
    * Your proposed solution to this problem, or your API design.
    * Alternative workarounds you have already considered.
3. Wait for Feedback: Let the project maintainers discuss the design with you. This should prevent you wasting time on code that
   may not align with the project's goals.

 ## Getting Started
 1. **Fork the repository** and clone it to your computer.
 2. **Install the Rust toolchain** using [rustup](https://rustup.rs).
 Ensure you have the latest version: `rustup update stable`
 3. Create a new branch for your work: `git checkout -b feature/my-feature-name`

## Development Workflow
Before submitting, please take the following steps:
* **Formatting:** Run `cargo fmt` to auto-format your code.
* **Linting:** Run `cargo clippy --all-targets --all-features`
* **Testing:** Run `cargo test`

## Documentation
1. If you add a new public API, structure or function:
   * Document Outer Doc Comments: Document items using /// for structures, functions, and modules.
   * Document Inner Doc Comments: Use //! at the very top of files to document crates or module-level overviews.
   * Leverage Common Headings: Organise long text using # Examples, #Errors, or # Panics sections.
   * Utilize Markdown Features: Write clear descriptions using standard markdown lists, paragraphs, and backticks (triple backticks) for code.
   * Apply Intradoc Links: Link to other types using [Type] format to leverage [Rusdoc's automated link generation].
   * Include runnable code examples in the documentation where applicable.
   * Commit your changes with clear, descriptive commit messages.
   * Write executable examples (doc tests): wrap code examples in standard markdown triple backticks. Ensure these code blocks are fully functional.
2. Update your branch to the latest `main` with `git pull origin main` to open a Pull Request (PR) against our main branch.
3. Link your PR to the approved feature request issue.
4. Push your changes to your fork and `[create a Pull Request on Github]`
5. Ensure your Pr description clearly states the problem solved


## Code Review Process
// TODO

## Recommended Best Practices
* API Guidelines: As a contributor you are encouraged to align with the Rust API Guidelines for idiomatic naming and design.

## Getting Help
// TODO

## Additional Resources
1. Google AI
2. Medium - iamprovidence
3. GitHub Docs: docs.github.com
4. Willow: willow voice.com: How to write a Good Pull Request Description: Developer's Guide for May 2026.
5. intersect-training.org: making Good PRs: Basic Pull Requests-(Last updated on 2026-06-22).
6. deployed.com: Pull Request Best Practices: A Complete Guide (2026).
7. viranaj.medium.com: Pull Request best Practices: A Comprehensive Guide for Developers.
8. Kara Luton: A Guide to Perfecting Pull Requests- DEV Community.







