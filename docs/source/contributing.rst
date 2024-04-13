Contributing Guidelines
=======================

Thank you for considering contributing to our project! We welcome contributions from everyone, and we appreciate your efforts to help improve our project. By contributing to this project, you agree to abide by our `Code of Conduct <CODE_OF_CONDUCT.md>`_.

Ways to Contribute
-------------------

There are many ways you can contribute to our project:

1. **Reporting Bugs:** If you encounter a bug or unexpected behavior, please `open an issue <https://github.com/Isy89/LBF/issues>`_ and provide detailed information about the problem.
2. **Suggesting Enhancements:** If you have ideas for new features or improvements, feel free to `open an issue <https://github.com/Isy89/LBF/issues>`_ to discuss them.
3. **Submitting Pull Requests:** If you'd like to contribute code changes, please fork our repository, create a new branch for your changes, and submit a pull request. Make sure to follow our `pull request template <PULL_REQUEST_TEMPLATE.md>`_ and adhere to our coding conventions.
4. **Improving Documentation:** You can help improve our project's documentation by fixing typos, clarifying explanations, or adding new content. Simply submit a pull request with your changes.
5. **Providing Feedback:** Your feedback is valuable to us! Whether it's about the project's usability, documentation, or any other aspect, feel free to share your thoughts by opening an issue or reaching out to us directly.

Getting Started
---------------

To start contributing to our project, follow these steps:

1. **Fork the Repository:** Click the "Fork" button on our repository's page to create your own copy of the project.
2. **Clone the Repository:** Clone your forked repository to your local machine using the following command::

   $ git clone https://github.com/Isy89/LBF.git

3. **Create a New Branch:** Create a new branch for your changes using a descriptive name. For example::

   $ git checkout -b feature/new-feature

4. **Make Changes:** Make your desired changes to the codebase, documentation, or any other aspect of the project.
5. **Test Your Changes:** Before submitting a pull request, make sure your changes work as expected and do not introduce any new issues. To To do this, run the tests.
Tests are run using pytest, which can be installed as follows::

    $ python -m pip install pytest
After intalling pytest, test can be run as follows::

    $ pytest -s tests
    
This should be done from within the LBF directory, which was created after cloning the git repo
6. **Bpdate the CHANGELOG:** update the CHANGELOG.md file explaining the changes or fixes.
7. **Bump the version:** bump the version in the setup.py file.
6. **Submit a Pull Request:** Push your branch to your forked repository and submit a pull request to our main repository.

Code of Conduct
---------------

Please note that by participating in this project, you are expected to uphold our `Code of Conduct <CODE_OF_CONDUCT.md>`. If you encounter any behavior that violates the code of conduct, please report it to us.

Enforcement Responsibilities
---------------------------

Community leaders are responsible for clarifying and enforcing our standards of acceptable behavior and will take appropriate and fair corrective action in response to any behavior that they deem inappropriate, threatening, offensive, or harmful.

Community leaders have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned with this Code of Conduct, and will communicate reasons for moderation decisions when appropriate.

Scope
-----

This Code of Conduct applies both within project spaces and in public spaces when an individual is representing the project or its community. Examples of representing a project or community include using an official project email address, posting via an official social media account, or acting as an appointed representative at an online or offline event.

Enforcement
-----------

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at `LBFextract@gmail.com`. All complaints will be reviewed and investigated and will result in a response that is deemed necessary and appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project's leadership.

Attribution
-----------

This Code of Conduct is adapted from the `Contributor Covenant <https://www.contributor-covenant.org>`_, version 2.1, available at `https://www.contributor-covenant.org/version/2/1/code_of_conduct.html`.

For answers to common questions about this code of conduct, see `FAQ <https://www.contributor-covenant.org/faq>`_.
