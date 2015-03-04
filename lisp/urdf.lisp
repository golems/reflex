;;;; -*- mode: lisp -*-
;;;;
;;;; Copyright (c) 2015, Rice University
;;;; All rights reserved.
;;;;
;;;; Author(s): Neil T. Dantam <ntd@rice.edu>
;;;;
;;;;   Redistribution and use in source and binary forms, with or
;;;;   without modification, are permitted provided that the following
;;;;   conditions are met:
;;;;   * Redistributions of source code must retain the above copyright
;;;;     notice, this list of conditions and the following disclaimer.
;;;;   * Redistributions in binary form must reproduce the above
;;;;     copyright notice, this list of conditions and the following
;;;;     disclaimer in the documentation and/or other materials provided
;;;;     with the distribution.
;;;;   * Neither the name of the copyright holder nor the names of its
;;;;     contributors may be used to endorse or promote products
;;;;     derived from this software without specific prior written
;;;;     permission.
;;;;   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
;;;;   CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
;;;;   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
;;;;   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
;;;;   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
;;;;   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
;;;;   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
;;;;   NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
;;;;   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
;;;;   HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
;;;;   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
;;;;   OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
;;;;   EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

(in-package :reflex)

;;; URDF RPY euler angles are fixed XYZ angles

;;; Translation of URDF to kinematic frames:
;;;  * Only joints matter
;;;  * When joints have rotated origins, add prior fixed frame for the origin.
;;;    - Ideally, the URDF source should be written to avoid this
;;;      since proper selection of axes avoids the need for these
;;;      intermediate frames

(defun slurp (file)
  (if (streamp file)
      (let ((string (make-string (file-length file))))
        (read-sequence string file)
        string)
      (with-open-file (s file :direction :input)
        (slurp s))))


(defstruct urdf-joint
  name
  rpy
  xyz
  type
  axis
  parent
  child)

(defun parse-string-vector (string)
  (when string
    (map 'list #'amino:parse-float (split-sequence:split-sequence #\Space string))))

(defun read-urdf (file)
  (let ((se (xmls:parse (slurp file))))
    (labels ((get-attr (name attrs)
               (second (assoc name attrs :test #'string=))))
    (destructuring-bind (top-elt top-attrs &rest sub-elts) se
      (let ((links)
            (joints))
        (dolist (x sub-elts)
          (destructuring-bind (elt attrs &rest sub) x
            (cond
              ((string= "link" elt)
               (push (get-attr "name" attrs)
                     links))
              ((string= "joint" elt)
               (let ((joint (make-urdf-joint
                             :name (get-attr "name" attrs)
                             :type (get-attr "type" attrs))))
                 (dolist (x sub)
                   (destructuring-bind (elt attrs &rest sub) x
                     (declare (ignore sub))
                     (cond
                       ((string= "origin" elt)
                        (setf (urdf-joint-rpy joint)
                              (parse-string-vector (get-attr "rpy" attrs)))
                        (setf (urdf-joint-xyz joint)
                              (parse-string-vector (get-attr "xyz" attrs))))
                       ((string= "parent" elt)
                        (setf (urdf-joint-parent joint) (get-attr "link" attrs)))
                       ((string= "child" elt)
                        (setf (urdf-joint-child joint) (get-attr "link" attrs)))
                       ((string= "axis" elt)
                        (setf (urdf-joint-axis joint)
                              (parse-string-vector (get-attr "xyz" attrs)))))))
                 (push joint joints))))))
        (values links joints))))))

(defun sort-frames (frames)
  (let* ((out)
         (hash (make-hash-table :test #'equal)))
    ;; initialize hash
    (loop for f in frames
       do (setf (gethash (frame-name f) hash)
                f))
    ;; process
    (labels (;(get-frame (name) (gethash name hash))
             (visit (parent)
               (push parent out)
               (loop with name = (frame-name parent)
                  for f in frames
                  when (string= (frame-parent f) name)
                  do (visit f))))
      ;; initialize stack with base tfs
      (loop for f in frames
         unless (frame-parent f)
         do (visit f)))
    (reverse out)))

(defun rpy->quaternion-list (rpy)
  (let ((q (amino:quaternion (amino:euler-zyx* (elt rpy 2)
                                                (elt rpy 1)
                                                (elt rpy 0)))))
    (list (amino:quaternion-x q)
          (amino:quaternion-y q)
          (amino:quaternion-z q)
          (amino:quaternion-w q))))

(defun urdf-joint-type-keyword (joint)
  (let ((s (string-downcase (urdf-joint-type joint))))
    (cond
      ((string= "revolute" s)
       :revolute)
      ((string= "prismatic" s)
       :revolute)
      ((string= "fixed" s)
       :fixed)
      (t
       (error "Unknown joint type ~A" s)))))

(defun parse-urdf (file)
  (multiple-value-bind (links joints)
      (read-urdf file)
    (declare (ignore links))
    (let ((link-parents (make-hash-table :test #'equal))
          (joint-parents (make-hash-table :test #'equal)))
      ;; find link parents
      (loop for j in joints
         for link = (urdf-joint-child j)
         for name = (urdf-joint-name j)
         when link
         do
           (setf (gethash link link-parents)
                 name))
      ;; find joint parents
      (loop for j in joints
         for link = (urdf-joint-parent j)
         for name = (urdf-joint-name j)
         for p = (gethash link link-parents)
         when p
         do (setf (gethash name joint-parents)
                  p))
      (let ((frames (loop for j in joints
                       for name = (urdf-joint-name j)
                       for parent = (gethash name joint-parents)
                       nconc
                         (cond
                           ;; fixed joint
                           ((string= (urdf-joint-type j) "fixed")
                            (list (make-frame name :fixed parent
                                              :quaternion (rpy->quaternion-list (urdf-joint-rpy j))
                                              :translation (urdf-joint-xyz j))))
                           ;; not rotated origin
                           ((equal (vector 0.0 0.0 0.0)
                                   (urdf-joint-rpy j))
                            (list (make-frame name (urdf-joint-type j) parent
                                              :axis (urdf-joint-axis j)
                                              :quaternion (rpy->quaternion-list (urdf-joint-rpy j))
                                              :translation (urdf-joint-xyz j))))
                           ;; rotated origin, need intermediate frame
                           (t
                            (let ((origin (concatenate 'string name "_origin")))
                              (list (make-frame origin :fixed parent
                                                :quaternion (rpy->quaternion-list (urdf-joint-rpy j))
                                                :translation (urdf-joint-xyz j))
                                    (make-frame name (urdf-joint-type j) origin
                                                :axis (urdf-joint-axis j)
                                                :translation '(0 0 0)
                                                :quaternion '(0 0 0 1)))))))))
        (sort-frames frames)))))
